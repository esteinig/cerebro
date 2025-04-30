use std::path::{Path, PathBuf};
use clap::{Args, ValueEnum};
use hf_hub::api::sync::Api;
use std::io::Write;

use tokenizers::Tokenizer;
use candle_examples::token_output_stream::TokenOutputStream;

use candle_core::{Device, Tensor};
use candle_core::quantized::gguf_file;

use candle_transformers::utils::apply_repeat_penalty;
use candle_transformers::models::quantized_llama as llama;
use candle_transformers::models::quantized_qwen2 as qwen2;
use candle_transformers::generation::{LogitsProcessor, Sampling};

use crate::error::GptError;

// Define a trait that abstracts over your two weight types
pub trait InferenceModel {
    /// run a single forward pass at a given position
    fn forward(&mut self, input: &Tensor, position: usize) -> Result<Tensor, GptError>;
}

impl InferenceModel for llama::ModelWeights {
    fn forward(&mut self, input: &Tensor, position: usize) -> Result<Tensor, GptError> {
        // delegate to the inherent method on ModelWeights
        Ok(llama::ModelWeights::forward(self, input, position)?)
    }
}

impl InferenceModel for qwen2::ModelWeights {
    fn forward(&mut self, input: &Tensor, position: usize) -> Result<Tensor, GptError> {
        // delegate to the inherent method on ModelWeights
        Ok(qwen2::ModelWeights::forward(self, input, position)?)
    }
}

pub struct TextGenerator {
    pub model: GeneratorModel,
    pub model_dir: PathBuf,
    pub force_download: bool
}

impl TextGenerator {
    pub fn new(model: GeneratorModel, model_dir: &Path, force_download: bool) -> Self {

        // Mixed precision GEMM kernels required for quantized Llama
        match model {
            GeneratorModel::DeepseekR1Llama8bQ4KM => {
                log::info!("Activate reduced precision GEMM kernels");
                candle_core::cuda::set_gemm_reduced_precision_f16(true);
                candle_core::cuda::set_gemm_reduced_precision_bf16(true);
            },
            _ => {},
        }

        Self {
            model,
            model_dir: model_dir.to_path_buf(),
            force_download
        }
    }
    pub fn run(
        &self,
        prompt: &str,
        raw_prompt: bool,
        sample_len: usize,
        split_prompt: bool,      // process prompt elements separately
        clean: bool,
        temperature: f64,
        seed: u64,               // seed for generating random samples
        top_k: Option<usize>,    // only sample among the top k samples
        top_p: Option<f64>,      // nucleus sampling probability cutoff
        repeat_penalty: f32,
        repeat_last_n: usize,
        log_info: bool
    ) -> Result<(), GptError> {

        let (model_path, tokenizer_path) = self.get_model()?;

        let device = candle_examples::device(false)?;
        let mut file = std::fs::File::open(&model_path)?;

        let start_reading_tensors = std::time::Instant::now();
        let gguf = gguf_file::Content::read(&mut file)
            .map_err(|e| e.with_path(&model_path))?;
        self.gguf_info(&gguf, &start_reading_tensors);

        let tokenizer = Tokenizer::from_file(tokenizer_path)?;
        let mut tos = TokenOutputStream::new(tokenizer);

        let prompt = if raw_prompt {
            prompt.to_string()
        } else {
            self.model.format_prompt(prompt)
        };

        if log_info {
            log::info!("Prompt is: \n\n{prompt}\n\n");
        }

        let tokens = tos
            .tokenizer()
            .encode(prompt, true)?;

        let tokens = tokens.get_ids();
        let to_sample = sample_len.saturating_sub(1);


        let mut logits_processor = self.build_logits_processor(
            temperature, 
            seed, 
            top_k, 
            top_p
        );

        match self.model {
            GeneratorModel::DeepseekR1Llama8bQ4KM => {

                let tokens = if tokens.len() + to_sample > llama::MAX_SEQ_LEN - 10 {
                    let to_remove = tokens.len() + to_sample + 10 - llama::MAX_SEQ_LEN;
                    tokens[tokens.len().saturating_sub(to_remove)..].to_vec()
                } else {
                    tokens.to_vec()
                };

                let mut model = llama::ModelWeights::from_gguf(
                    gguf, 
                    &mut file, 
                    &device
                )?;

                self.generate(
                    &mut model,
                    &device,
                    &tokens,
                    &mut tos,
                    &mut logits_processor,
                    sample_len,
                    repeat_penalty,
                    repeat_last_n,
                    split_prompt,
                    log_info,
                    clean
                )?;
            }
            GeneratorModel::DeepseekR1Qwen7bQ4KM 
                | GeneratorModel::DeepseekR1Qwen7bQ80
                | GeneratorModel::DeepseekR1Qwen7bF16
                | GeneratorModel::DeepseekR1Qwen14bQ2KL
                | GeneratorModel::DeepseekR1Qwen14bQ4KM
                | GeneratorModel::DeepseekR1Qwen14bQ80
                | GeneratorModel::DeepseekR1Qwen14bF16 
             => {

                let mut model = qwen2::ModelWeights::from_gguf(
                    gguf, 
                    &mut file, 
                    &device
                )?;

                self.generate(
                    &mut model,
                    &device,
                    tokens,
                    &mut tos,
                    &mut logits_processor,
                    to_sample,
                    repeat_penalty,
                    repeat_last_n,
                    split_prompt,
                    log_info,
                    clean
                )?;
            }
        };


        Ok(())
        

    }
    pub fn generate(
        &self,
        model: &mut dyn InferenceModel,
        device: &Device,
        tokens: &[u32],
        tos: &mut TokenOutputStream,
        logits_processor: &mut LogitsProcessor,
        to_sample: usize,
        repeat_penalty: f32,
        repeat_last_n: usize,
        split_prompt: bool,
        log_info: bool,
        clean: bool
    ) -> Result<(), GptError> {

        // Holds all tokens generated
        let mut all_tokens = vec![];
        
        // Process the prompt input
        let (first_token, prompt_dt) = self.process_prompt(
            &mut *model,
            &tokens,
            split_prompt,
            &device,
            logits_processor,
        )?;

        // Hold the token and print it out
        all_tokens.push(first_token);
        if let Some(t) = tos.next_token(first_token)? {
            print!("{t}");
            std::io::stdout().flush()?;
        }

        let eos = self.model.get_eos_token(tos)?;

        // Process the main sample loop
        let (sampled, gen_dt) = self.generate_qwen(
            &mut *model,
            first_token,
            &mut all_tokens,
            logits_processor,
            to_sample,
            tokens.len(),
            repeat_penalty,
            repeat_last_n,
            eos,
            tos,
            &device,
        )?;

        // Flush any trailing subwords and final newline
        if let Some(rest) = tos.decode_rest().map_err(candle_core::Error::msg)? {
            print!("{rest}");
        }
        std::io::stdout().flush()?;

        println!("\n");

        if log_info {
            log::info!(
                "{:4} prompt tokens processed @ {:.2} token/s",
                tokens.len(),
                tokens.len() as f64 / prompt_dt.as_secs_f64(),
            );
            log::info!(
                "{:4} tokens generated @ {:.2} token/s",
                sampled,
                sampled as f64 / gen_dt.as_secs_f64(),
            );
        }
        
        Ok(())
    }
    /// Run the post-prompt generation loop for Qwen-based models
    pub fn generate_qwen(
        &self,
        model: &mut dyn InferenceModel,
        mut next_token: u32,
        all_tokens: &mut Vec<u32>,
        logits_processor: &mut LogitsProcessor,
        to_sample: usize,
        prompt_len: usize,
        repeat_penalty: f32,
        repeat_last_n: usize,
        eos_token: u32,
        tos: &mut TokenOutputStream,
        device: &Device,
    ) -> Result<(usize, std::time::Duration), GptError> {

        let start = std::time::Instant::now();
        let mut sampled = 0;

        for i in 0..to_sample {

            // Forward one token
            let input = Tensor::new(&[next_token], device)?.unsqueeze(0)?;
            let mut logits = model.forward(&input, prompt_len + i)?.squeeze(0)?;

            // Optional repeat-penalty
            if repeat_penalty != 1.0 {
                let begin = all_tokens.len().saturating_sub(repeat_last_n);
                logits = apply_repeat_penalty(&logits, repeat_penalty, &all_tokens[begin..])?;
            }

            // Sample, record, print & break on EOS
            next_token = logits_processor.sample(&logits)?;
            all_tokens.push(next_token);

            if let Some(text) = tos.next_token(next_token)? {
                print!("{}", text);
                std::io::stdout().flush()?;
            }

            sampled += 1;
            if next_token == eos_token {
                break;
            }
        }

        Ok((sampled, start.elapsed()))
    }
    pub fn process_prompt(
        &self,
        model: &mut dyn InferenceModel,
        tokens: &[u32],
        split_prompt: bool,
        device: &Device,
        logits_processor: &mut LogitsProcessor,
    ) -> Result<(u32, std::time::Duration), GptError> {
        let start = std::time::Instant::now();
        
        let tok = if !split_prompt {
            let input = Tensor::new(tokens, device)?.unsqueeze(0)?;
            let logits = model.forward(&input, 0)?.squeeze(0)?;
            logits_processor.sample(&logits)?
        } else {
            let mut last = 0;
            for (pos, &tok) in tokens.iter().enumerate() {
                let single = Tensor::new(&[tok], device)?.unsqueeze(0)?;
                let logits = model.forward(&single, pos)?.squeeze(0)?;
                last = logits_processor.sample(&logits)?;
            }
            last
        };

        Ok((tok, start.elapsed()))
    }
    pub fn build_logits_processor(&self, temperature: f64, seed: u64, top_k: Option<usize>, top_p: Option<f64>) -> LogitsProcessor {

        let sampling = if temperature <= 0. {
            Sampling::ArgMax
        } else {
            match (top_k, top_p) {
                (None, None) => Sampling::All { temperature },
                (Some(k), None) => Sampling::TopK { k, temperature },
                (None, Some(p)) => Sampling::TopP { p, temperature },
                (Some(k), Some(p)) => Sampling::TopKThenTopP { k, p, temperature },
            }
        };
        LogitsProcessor::from_sampling(seed, sampling)
    }
    pub fn gguf_info(&self, gguf: &gguf_file::Content, start: &std::time::Instant) {

        let mut total_size_in_bytes = 0;
        for (_, tensor) in gguf.tensor_infos.iter() {
            let elem_count = tensor.shape.elem_count();
            total_size_in_bytes +=
                elem_count * tensor.ggml_dtype.type_size() / tensor.ggml_dtype.block_size();
        }
        log::info!(
            "loaded {:?} tensors ({}) in {:.2}s",
            gguf.tensor_infos.len(),
            &format_size(total_size_in_bytes),
            start.elapsed().as_secs_f32(),
        );
    }
    pub fn get_model(&self) -> Result<(PathBuf, PathBuf), GptError> {

        let model_path = self.model_dir.join(&self.model.model_file());
        let tokenizer_path = self.model_dir.join(&self.model.tokenizer_file());
        
        let model_file = if model_path.exists() & !self.force_download {
            model_path.clone()
        } else {
            self.model.save_model(&self.model_dir)?
        };

        let tokenizer_file = if tokenizer_path.exists() & !self.force_download {
            tokenizer_path.clone()
        } else {
            self.model.save_tokenizer(&self.model_dir)?
        };

        Ok((model_file, tokenizer_file))
    }
}

#[derive(Clone, Debug, Copy, PartialEq, Eq, ValueEnum)]
pub enum GeneratorModel {
    #[value(name = "deepseekr1-llama8b-q4-km")]
    DeepseekR1Llama8bQ4KM,
    
    #[value(name = "deepseekr1-qwen7b-q4-km")]
    DeepseekR1Qwen7bQ4KM,
    #[value(name = "deepseekr1-qwen7b-q8-0")]
    DeepseekR1Qwen7bQ80,
    #[value(name = "deepseekr1-qwen7b-f16")]
    DeepseekR1Qwen7bF16,

    #[value(name = "deepseekr1-qwen14b-q2-kl")]
    DeepseekR1Qwen14bQ2KL,
    #[value(name = "deepseekr1-qwen14b-q4-km")]
    DeepseekR1Qwen14bQ4KM,
    #[value(name = "deepseekr1-qwen14b-q8-0")]
    DeepseekR1Qwen14bQ80,
    #[value(name = "deepseekr1-qwen14b-f16")]
    DeepseekR1Qwen14bF16,
}

impl GeneratorModel {

    /// Download and save the GGUF model file as `{model_name}.{ext}`
    pub fn save_model(&self, outdir: &Path) -> Result<PathBuf, GptError> {
        let src = self.download_model()?;
        std::fs::create_dir_all(outdir)?;
        let dest = outdir.join(&self.model_file());
        std::fs::copy(&src, &dest)?;
        Ok(dest)
    }
    /// Download and save the tokenizer to `{model_name}.tokenizer.json`
    pub fn save_tokenizer(&self, outdir: &Path) -> Result<PathBuf, GptError> {
        let src = self.download_tokenizer()?;
        std::fs::create_dir_all(outdir)?;
        let dest = outdir.join(&self.tokenizer_file());
        std::fs::copy(&src, &dest)?;
        Ok(dest)
    }
    pub fn download_tokenizer(&self) -> Result<PathBuf, GptError> {
        log::info!("Downloading tokenizer...");
        let api = Api::new()?;
        let repo = self.tokenizer_repository();
        let api = api.model(repo.to_string());
        let tokenizer_path = api.get("tokenizer.json")?;
        Ok(tokenizer_path)
    }
    pub fn download_model(&self) -> Result<PathBuf, GptError> {
        log::info!("Downloading model weights...");
        let api = Api::new()?;
        let repo = hf_hub::Repo::with_revision(
            self.model_repository().to_string(),
            hf_hub::RepoType::Model,
            self.model_revision().to_string(),
        );
        let model_path = api.repo(repo).get(
            self.model_config()
        )?;
        Ok(model_path)
    }
    pub fn tokenizer_repository(&self) -> &'static str {
        match self {
            GeneratorModel::DeepseekR1Llama8bQ4KM 
                => "deepseek-ai/DeepSeek-R1-Distill-Llama-8B",
            GeneratorModel::DeepseekR1Qwen7bQ4KM  
                | GeneratorModel::DeepseekR1Qwen7bQ80
                | GeneratorModel::DeepseekR1Qwen7bF16  
                => "deepseek-ai/DeepSeek-R1-Distill-Qwen-7B",
            GeneratorModel::DeepseekR1Qwen14bQ2KL
                | GeneratorModel::DeepseekR1Qwen14bQ4KM
                | GeneratorModel::DeepseekR1Qwen14bQ80
                | GeneratorModel::DeepseekR1Qwen14bF16
                => "deepseek-ai/DeepSeek-R1-Distill-Qwen-14B",
        }
    }
    pub fn model_repository(&self) -> &'static str {
        match self {
            GeneratorModel::DeepseekR1Llama8bQ4KM 
                => "unsloth/DeepSeek-R1-Distill-Llama-8B-GGUF",
            GeneratorModel::DeepseekR1Qwen7bQ4KM 
                | GeneratorModel::DeepseekR1Qwen7bQ80
                | GeneratorModel::DeepseekR1Qwen7bF16
                => "unsloth/DeepSeek-R1-Distill-Qwen-7B-GGUF",
            GeneratorModel::DeepseekR1Qwen14bQ2KL
                | GeneratorModel::DeepseekR1Qwen14bQ4KM
                | GeneratorModel::DeepseekR1Qwen14bQ80
                | GeneratorModel::DeepseekR1Qwen14bF16
                => "unsloth/DeepSeek-R1-Distill-Qwen-14B-GGUF",
        }
    }
    pub fn model_config(&self) -> &'static str {
        match self {
            GeneratorModel::DeepseekR1Llama8bQ4KM    => "DeepSeek-R1-Distill-Llama-8B-Q4_K_M.gguf",
            GeneratorModel::DeepseekR1Qwen7bQ4KM => "DeepSeek-R1-Distill-Qwen-7B-Q4_K_M.gguf",
            GeneratorModel::DeepseekR1Qwen7bQ80 => "DeepSeek-R1-Distill-Qwen-7B-Q8_0.gguf",
            GeneratorModel::DeepseekR1Qwen7bF16 => "DeepSeek-R1-Distill-Qwen-7B-F16.gguf",
            GeneratorModel::DeepseekR1Qwen14bQ2KL => "DeepSeek-R1-Distill-Qwen-14B-Q2_K_L.gguf",
            GeneratorModel::DeepseekR1Qwen14bQ4KM => "DeepSeek-R1-Distill-Qwen-14B-Q4_K_M.gguf",
            GeneratorModel::DeepseekR1Qwen14bQ80 => "DeepSeek-R1-Distill-Qwen-14B-Q8_0.gguf",
            GeneratorModel::DeepseekR1Qwen14bF16 => "DeepSeek-R1-Distill-Qwen-14B-F16.gguf",
        }
    }
    pub fn model_revision(&self) -> &'static str {
        "main"
    }
    // For writing files to disk
    pub fn model_name(&self) -> &'static str {
        match self {
            GeneratorModel::DeepseekR1Llama8bQ4KM    => "deepseekr1-llama8b-q4-km",
            GeneratorModel::DeepseekR1Qwen7bQ4KM => "deepseekr1-qwen7b-q4-km",
            GeneratorModel::DeepseekR1Qwen7bQ80 => "deepseekr1-qwen7b-q8-0",
            GeneratorModel::DeepseekR1Qwen7bF16 => "deepseekr1-qwen7b-f16",
            GeneratorModel::DeepseekR1Qwen14bQ2KL => "deepseekr1-qwen14b-q2-kl",
            GeneratorModel::DeepseekR1Qwen14bQ4KM => "deepseekr1-qwen14b-q4-km",
            GeneratorModel::DeepseekR1Qwen14bQ80 => "deepseekr1-qwen14b-q8-0",
            GeneratorModel::DeepseekR1Qwen14bF16 => "deepseekr1-qwen14b-f16",
        }
    }
    pub fn model_file(&self) -> PathBuf {
        let model_config = PathBuf::from(
            self.model_config()
        );
        let ext = model_config
            .extension()
            .and_then(|e| e.to_str())
            .unwrap_or("config");
        PathBuf::from(format!("{}.{}", self.model_name(), ext))
    }

    pub fn tokenizer_file(&self) -> PathBuf {
        PathBuf::from(format!("{}.tokenizer.json", self.model_name()))
    }

    pub fn get_eos_token(&self, 
        tos: &mut TokenOutputStream,
    ) -> Result<u32, GptError> {

        // Get the  end of sentence token
        let eos_token = match self {
            GeneratorModel::DeepseekR1Llama8bQ4KM 
                | GeneratorModel::DeepseekR1Qwen7bQ4KM 
                | GeneratorModel::DeepseekR1Qwen7bQ80
                | GeneratorModel::DeepseekR1Qwen7bF16
                | GeneratorModel::DeepseekR1Qwen14bQ2KL
                | GeneratorModel::DeepseekR1Qwen14bQ4KM
                | GeneratorModel::DeepseekR1Qwen14bQ80
                | GeneratorModel::DeepseekR1Qwen14bF16  => "<｜end▁of▁sentence｜>"
        };
        let eos = *tos.tokenizer()
            .get_vocab(true)
            .get(eos_token)
            .ok_or(GptError::EosTokenNotInVocabulary(eos_token.to_string()))?;

        Ok(eos)
    }

    pub fn is_deepseek_qwen(&self) -> bool {
        match self {
            GeneratorModel::DeepseekR1Llama8bQ4KM => false,
            GeneratorModel::DeepseekR1Qwen7bQ4KM 
            | GeneratorModel::DeepseekR1Qwen7bQ80
            | GeneratorModel::DeepseekR1Qwen7bF16
            | GeneratorModel::DeepseekR1Qwen14bQ2KL
            | GeneratorModel::DeepseekR1Qwen14bQ4KM
            | GeneratorModel::DeepseekR1Qwen14bQ80
            | GeneratorModel::DeepseekR1Qwen14bF16 => true,
        }
    }
    pub fn is_deepseek_llama(&self) -> bool {
        match self {
            GeneratorModel::DeepseekR1Llama8bQ4KM => true,
            GeneratorModel::DeepseekR1Qwen7bQ4KM 
            | GeneratorModel::DeepseekR1Qwen7bQ80
            | GeneratorModel::DeepseekR1Qwen7bF16
            | GeneratorModel::DeepseekR1Qwen14bQ2KL
            | GeneratorModel::DeepseekR1Qwen14bQ4KM
            | GeneratorModel::DeepseekR1Qwen14bQ80
            | GeneratorModel::DeepseekR1Qwen14bF16 => false
        }
    }

    pub fn format_prompt(&self, prompt: &str) -> String {
        if self.is_deepseek_qwen() {
            format!("<｜User｜>{prompt}<｜Assistant｜>")
        } else if self.is_deepseek_llama() {
            format!("<｜user｜>{prompt}<｜assistant｜>")  // llama distillation only works with non-capitalized tags?
        } else {
            prompt.to_string()
        }
    }
}


fn format_size(size_in_bytes: usize) -> String {
    if size_in_bytes < 1_000 {
        format!("{}B", size_in_bytes)
    } else if size_in_bytes < 1_000_000 {
        format!("{:.2}KB", size_in_bytes as f64 / 1e3)
    } else if size_in_bytes < 1_000_000_000 {
        format!("{:.2}MB", size_in_bytes as f64 / 1e6)
    } else {
        format!("{:.2}GB", size_in_bytes as f64 / 1e9)
    }
}


#[derive(Args, Debug, Clone)]
pub struct TextGeneratorArgs {

    /// Text generation model.
    #[arg(long, default_value = "deepseekr1-qwen7b-q4-km")]
    pub model: GeneratorModel,

    /// Input user prompt.
    #[arg(long)]
    pub prompt: String,

    /// Model file and download directory.
    #[arg(long, default_value=".")]
    pub model_dir: PathBuf,

    /// Force download of model and tokenizer files if they do not exist in the model directory.
    #[arg(long)]
    pub force_download: bool,

    /// Force raw prompt instead of adding model tokens to input prompt.
    #[arg(long)]
    pub raw_prompt: bool,

    /// The length of the sample to generate (in tokens).
    #[arg(short = 'n', long, default_value_t = 1000)]
    pub sample_len: usize,

    /// The temperature used to generate samples, use 0 for greedy sampling.
    #[arg(long, default_value_t = 0.8)]
    pub temperature: f64,

    /// Nucleus sampling probability cutoff.
    #[arg(long)]
    pub top_p: Option<f64>,

    /// Only sample among the top K samples.
    #[arg(long)]
    pub top_k: Option<usize>,

    /// The seed to use when generating random samples.
    #[arg(long, default_value_t = 299792458)]
    pub seed: u64,

    /// Process prompt elements separately.
    #[arg(long)]
    pub split_prompt: bool,

    /// Clean output prompt e.g. by removing <think> tags for Deepseek.
    #[arg(long, short='c')]
    pub clean: bool,

    /// Penalty to be applied for repeating tokens, 1. means no penalty.
    #[arg(long, default_value_t = 1.1)]
    pub repeat_penalty: f32,

    /// The context size to consider for the repeat penalty.
    #[arg(long, default_value_t = 64)]
    pub repeat_last_n: usize,

    /// Log additional information.
    #[arg(long)]
    pub log_info: bool,

}