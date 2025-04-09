// main.rs

use anyhow::Result;
use cerebro_client::client::CerebroClient;
use cerebro_model::api::cerebro::schema::{CerebroIdentifierSchema, MetaGpConfig};
use cerebro_pipeline::taxa::filter::TaxonFilterConfig;
use cerebro_pipeline::taxa::taxon::Taxon;
use petgraph::graph::{Graph, NodeIndex};
use petgraph::visit::EdgeRef;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::ops::Deref;
use std::path::{Path, PathBuf};
use tokio::process::Command;
use async_openai::{config::OpenAIConfig, types::*};
use colored::Colorize;
use crate::error::GptError;


//
// === Refined Question Types ===
//

#[derive(Debug, Deserialize)]
#[serde(untagged)]
pub enum Question {
    Simple(String),
    Detailed {
        prompt: String,
        context_description: Option<String>,
        instructions: Option<String>,
    },
}

impl Question {
    pub fn to_prompt(&self) -> String {
        match self {
            Question::Simple(s) => s.clone(),
            Question::Detailed {
                prompt,
                context_description,
                instructions,
            } => {
                let mut combined = String::new();
                combined.push_str(&format!("{}\n", prompt));
                if let Some(context) = context_description {
                    combined.push_str(&format!("Context: {}\n", context));
                }
                if let Some(instr) = instructions {
                    combined.push_str(&format!("Instructions: {}\n", instr));
                }
                combined
            }
        }
    }
}

#[derive(Clone, Debug, Deserialize, Serialize, clap::ValueEnum)]
pub enum GptStrategy {
    OneShotResponse,
    DecisionTreeChat,
    DecisionTreeResponse,
}

pub struct ThresholdCandidates {
    pub above_threshold: Option<Vec<Taxon>>,
    pub below_threshold: Option<Vec<Taxon>>,
    pub target_list: Option<Vec<Taxon>>
}

impl ThresholdCandidates {
    pub fn from_above_threshold(above_threshold: Vec<Taxon>) -> Self {
        Self { above_threshold: Some(above_threshold), below_threshold: None, target_list: None }
    }
    pub fn from_below_threshold(below_threshold: Vec<Taxon>) -> Self {
        Self { above_threshold: None, below_threshold: Some(below_threshold), target_list: None }
    }
    pub fn from_target_list(target_list: Vec<Taxon>) -> Self {
        Self { above_threshold: None, below_threshold: None, target_list: Some(target_list) }
    }
    /// Returns a formatted string summary of candidate taxa for each threshold.
    /// - If a candidate group is `None`, it is omitted.
    /// - If a candidate group is `Some` but empty, "no taxa called" is shown.
    /// - Otherwise, the header is printed along with the list of taxon names.
    pub fn to_str(&self) -> String {
        let mut output = String::new();

        // Process above threshold candidates
        if let Some(ref above) = self.above_threshold {
            output.push_str("Above threshold candidate taxa: ");
            if above.is_empty() {
                output.push_str("no taxa called\n\n");
            } else {
                // Join the names of the taxa with commas.
                let taxa: Vec<String> = above.iter().map(|taxon| format!("{taxon}")).collect();
                output.push_str(&taxa.join("\n"));
            }
            output.push('\n');
        }

        // Process below threshold candidates
        if let Some(ref below) = self.below_threshold {
            output.push_str("Below threshold candidate taxa: ");
            if below.is_empty() {
                output.push_str("no taxa called\n\n");
            } else {
                let taxa: Vec<String> = below.iter().map(|taxon| format!("{taxon}")).collect();
                output.push_str(&taxa.join("\n"));
            }
        }

        // Process target list candidate taxa
        if let Some(ref target) = self.target_list {
            output.push_str("Target threshold candidate taxa: ");
            if target.is_empty() {
                output.push_str("no taxa called\n\n");
            } else {
                let taxa: Vec<String> = target.iter().map(|taxon| format!("{taxon}")).collect();
                output.push_str(&taxa.join("\n"));
            }
        }

        output
    }
}


pub enum DataBackground {
    CerebroFilter
}
impl DataBackground {
    pub fn get_default(&self) -> String {
        match self {
            DataBackground::CerebroFilter => String::from("
                We conducted metagenomic sequencing for pathogen detection and diagnosis (Illumina PE, RNA and DNA libraries on NextSeq) - a 'needle in a haystack' problem. Filtering the taxonomic profiling data from 
                the Cerebro pipeline produced three subsets of the same dataset: above threshold (high-confidence, high abundance, specific but less sensitive), below threshold (sensitive but less specific, often 
                contamination and low-level pathogen abundance), and a target pathogen list of high priority pathogens of interest, in this case vertebrate viruses (very sensitive but much less specific, 
                often misalignments reported from viral reference queries low level,  filter on alkignment and k-mer evidence required imposed nevertheless). Cerebro in general uses multiple profiling methods 
                (alignment, k-mer, assembly) for pathogen detection. Cerebro uses k-mer classifiers (Kraken2+Bracken, Metabuli, Ganon2), read alignment (Vircov with bowtie2 against ICTV viral references only) 
                and metagenome assembly (Megahit) and BLAST queries of contigs - all methods use the same reference database Cipher which consists of ICTV, GTDB, EuPath, WormBase, FungiDB, etc,
                as well as the human reference genome (CHM13v2). There is a risk of real contamination from sampling, processing, the lab environment, as well as from reference database genome contamination 
                in particular in eukaryotic genomes as well as variable peformance sensitivties/specificities of the classifiers and aligners involved in making a taxon call. Host aneuploidy detection is conducted 
                in some sample DNA libraries which involves detection of abnormal copy number variation (CNV) of large chromosomal segments across the host genome as tumors and cancer can be a differential diagnosis
                for some infectious symptoms, and the metagenomics assay sequences a high amount of host background nucleic acid.
            "),
        }
    } 
}

#[derive(Clone, Debug, Deserialize, Serialize, clap::ValueEnum)]
pub enum ClinicalContext {
    Csf,
    Eye,
}

impl ClinicalContext {
    pub fn text(&self) -> String {
        match self {
            ClinicalContext::Csf => String::from("
                CSF sample from a patient with neuroinflammatory or neurological symptoms. Rare cases of infectious agents common in this sample type and clinical presentation should be considered if detected. Consider skin microbiome contamination from sampling site and other sources of contamination from handling of sterile samples.
            "),
            ClinicalContext::Eye => String::from("
                Vitreous fluid sample from a patient with ocular infection or neurological symptoms. Unusual bacterial species at high to medium abundance should be considered if detected. Consider skin microbiome contamination from sampling site and other sources of contamination from handling of sterile samples.
            "),
        }
    }
    pub fn with_default(&self, header: String, text: String) -> String {
        match self {
            ClinicalContext::Csf => format!("{}\n\n{header}==={text}", self.text()),
            ClinicalContext::Eye => format!("{}\n\n{header}==={text}", self.text()),
        }
    }
}

//
// === Decision Tree Definitions ===
//

#[derive(Debug, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum CheckType {
    LlmEval,
    LlmDiagnosticEval,
    LlmJsonEval,
    AboveThresholdQuery,
    BelowThresholdQuery,
    TargetThresholdQuery,
    AneuploidyQuery,
    TaxaHistoryQuery,
    FlightCheckNonInfectious
}

#[derive(Debug, Deserialize)]
pub struct TreeNode {
    question: Option<Question>,
    check: Option<CheckType>,
    true_node: Option<String>,
    false_node: Option<String>,
    next: Option<String>,
    final_node: Option<bool>,
}

pub type DecisionTree = HashMap<String, TreeNode>;

//
// === Diagnostic History Structures ===
//

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct DiagnosticStep {
    pub node_key: String,
    pub question: Option<String>,
    pub answer: Option<String>,
}

//
// === Evaluation Details Structure ===
//

#[derive(Debug, Serialize, Deserialize)]
pub struct EvaluationDetails {
    pub taxa_contamination: Vec<String>,
    pub taxa_pathogens: Vec<String>,
    pub contamination_reason: String,
    pub pathogen_reason: String,
    pub positive_infection: Option<bool>,
}


//
// === Final Diagnostic Result Structure ===
//

#[derive(Debug, Serialize, Deserialize, PartialEq)]
pub enum Diagnosis {
    Infectious,
    InfectiousReview,
    NonInfectious,
    NonInfectiousReview,
    Tumor,
    Unknown
}

impl std::fmt::Display for Diagnosis {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let s = match self {
            Diagnosis::Infectious => "infectious",
            Diagnosis::InfectiousReview => "infectious-review",
            Diagnosis::NonInfectious => "non-infectious",
            Diagnosis::NonInfectiousReview => "non-infectious-review",
            Diagnosis::Tumor => "tumor",
            Diagnosis::Unknown => "unknown",
        };
        write!(f, "{}", s)
    }
}


#[derive(Debug, Serialize, Deserialize)]
pub struct DiagnosticResult {
    pub diagnosis: Diagnosis,
    pub candidate: Option<String>,
    pub candidate_description: Option<String>,
    pub pathogen: Option<String>,
    pub pathogen_description: Option<String>,
    pub reason_non_infectious: Option<String>,
    pub clinical_context: String,
    pub diagnostic_history_shorthand: Vec<String>,
    pub diagnostic_history_annotated: Vec<DiagnosticStep>,
    pub taxa: ThresholdTaxa
}
impl DiagnosticResult {
    pub fn to_json(&self, path: &Path) -> Result<(), GptError> {
        let agent_state = serde_json::to_string_pretty(self).map_err(|err| GptError::SerdeJsonError(err))?;
        let mut writer = BufWriter::new(File::create(path)?);
        write!(writer, "{agent_state}")?;
        Ok(())
    }
    pub fn from_json<P: AsRef<Path>>(path: P) -> Result<Self, GptError> {
        let data = std::fs::read_to_string(path)?;
        let result = serde_json::from_str::<DiagnosticResult>(&data)?;
        Ok(result)
    }
}

//
// === Agent State (Memory) ===
//

#[derive(Debug, Serialize, Deserialize)]
pub struct AgentState {
    pub memory: HashMap<String, String>,
    pub history: Vec<DiagnosticStep>,
    pub shorthand: Vec<String>,
}

impl AgentState {
    fn new() -> Self {
        AgentState {
            memory: HashMap::new(),
            history: Vec::new(),
            shorthand: Vec::new(),
        }
    }

    fn log(&mut self, key: &str, value: &str) {
        self.memory.insert(key.to_string(), value.to_string());
    }

    fn add_history(&mut self, node_key: &str, question: Option<String>, answer: Option<&str>) {
        self.shorthand.push(node_key.to_string());
        self.history.push(DiagnosticStep {
            node_key: node_key.to_string(),
            question: question.map(|s| s.to_string()),
            answer: answer.map(|s| s.to_string()),
        });
    }
    fn get_history(&mut self, node_key: &str) -> Option<&DiagnosticStep> {
        self.history.iter().find(|s| s.node_key == node_key)
    }

    pub fn to_json(&mut self, path: &Path) -> Result<(), GptError> {
        let agent_state = serde_json::to_string_pretty(self).map_err(|err| GptError::SerdeJsonError(err))?;
        let mut writer = BufWriter::new(File::create(path)?);
        write!(writer, "{agent_state}")?;
        Ok(())
    }
}

//
// === Diagnostic Agent ===
//

#[derive(Debug, Serialize, Deserialize)]
pub struct ThresholdTaxa {
    pub above: Vec<Taxon>,
    pub below: Vec<Taxon>,
    pub target: Vec<Taxon>
}
impl ThresholdTaxa {
    pub fn from_agent_memory(memory: &HashMap<String, String>) -> Result<Self, GptError> {

        let above: Vec<Taxon> = serde_json::from_str(
            memory.get("data__AboveThresholdQuery__CandidateTaxa")
                  .map_or("[]", |r| r.deref())
        )?;

        let below: Vec<Taxon> = serde_json::from_str(
            memory.get("data__BelowThresholdQuery__CandidateTaxa")
                  .map_or("[]", |r| r.deref())
        )?;


        let target: Vec<Taxon> = serde_json::from_str(
            memory.get("data__TargetThresholdQuery__CandidateTaxa")
                  .map_or("[]", |r| r.deref())
        )?;

        Ok(Self {
            above,
            below,
            target
        })
    }
    pub fn to_json(&self, path: &Path) -> Result<(), GptError> {
        let mut writer = BufWriter::new(File::create(path)?);
        write!(&mut writer, "{}", serde_json::to_string_pretty(self)?)?;
        Ok(())
    }
}

pub struct DiagnosticAgent {
    pub tree: DecisionTree,
    pub cerebro_client: CerebroClient,
    pub llm_client: async_openai::Client<OpenAIConfig>,
    pub state: AgentState,
    pub knowledge_graph: Graph<String, String>,
    pub model: String,
    pub diagnostic_memory: bool,
    pub contam_history: bool
}

impl DiagnosticAgent {
    pub async fn new(cerebro_client: CerebroClient, model: String, diagnostic_memory: bool, contam_history: bool) -> Result<Self> {
        
        let decision_tree_json = r#"
        {
            "start": {
                "question": {
                    "prompt": "Does the segmental aneuploidy analysis (if conducted) indicate host tumor DNA in the sample? Copy number variation (CNV) along large segments of chromosomes are indicated in the scatter plot analysis and data extraction. Large variation (a braodly scattered) usually indicates low human genome coverage.",
                    "instructions": "Answer 'yes' if the data indicates a tumor signal, 'no' if a tumor-signal is not supported or the analysis was not conducted."
                },
                "check": "aneuploidy_query",
                "true_node": "diagnose_tumor",
                "false_node": "check_above_threshold"
            },
            "check_above_threshold": {
                "question": {
                    "prompt": "Based on the evidence synthesis from above-threshold metagenomic profiling, make a diagnosis of infectious or non-infectious in the context of the metagenomics assay, sources of contamination, microbial profile, patient clinical history, and sample type. If selection is not clear from the evidence synthesis or data is not avilable for this filter, you should select 'non-infectious'. Consider the presence of contaminant taxa at above threshold levels particularly from sample sites or lab environment. If pathogens known to commonly infect humans are detected you should select 'infectious'. If viral pathogens were detected from the target list and a detected virus is known to infect humans, consider this sample 'infectious'.",
                    "instructions": "Answer 'yes' if the above-threshold data supports infection, 'no' if it does not or if no taxa were called."
                },
                "check": "above_threshold_query",
                "true_node": "diagnose_infectious",
                "false_node": "check_below_threshold"
            },
            "check_below_threshold": {
                "question": {
                    "prompt": "Based on the evidence synthesis from below-threshold metagenomic profiling, make a diagnosis of infectious or non-infectious in the context of the metagenomics assay, sources of contamination, microbial profile, patient clinical history, sample type. If selection is not clear from the evidence synthesis or data is not available for this filter, you should select 'non-infectious'. If pathogens known to infect humans are detected and a detected target taxon has reasonable evidence you should select 'infectious'. This section sometimes contains low-level bacterial assemblages - if you are making this call based on multiple pathogen candidates, consider the overall diversity and alternative reasoning around contamination. If a taxon has both a pathogen and contaminant role -- e.g. a pathogen that is a commensal in sample or sample adjacent sites of the human body, deriving for example from skin or lab environment sequenced in the metagenomics assay -- deprioritize the taxon as a pathogen candidate especially if itoccurss at low abundance in the sample and among the microbial profile called in this section.",
                    "instructions": "Answer 'yes' if the below-threshold data supports infectious diagnosis; 'no' if it does not or if no taxa were called."
                },
                "check": "below_threshold_query",
                "next": "check_target_threshold"
            },
            "check_target_threshold": {
                "question": {
                    "prompt": "Based on the synthesis of target-list metagenomic profiling evidence, make a diagnosis of infectious or non-infectious in the context of the metagenomics assay, sources of contamination, microbial profile, patient clinical history, and sample type. If selection is not clear from the evidence synthesis or data is not avilable for this filter, you should select 'non-infectious'. If viral targets were detected from the target list and a detected target is known to infect humans, consider this sample 'infectious'. Other taxa from prokaryotic and eukaryotic domains should be considered with caution and only selected with reasonably strong evidence from the taxon calling in relation to clinical context of the patient from whom this sample derives.",
                    "instructions": "Answer 'yes' if the target list supports infectious diagnosis; 'no' if it does not or if no taxa were called."
                },
                "check": "target_threshold_query",
                "next": "integrate_below_target_evidence"
            },
            "integrate_below_target_evidence": {
                "question": {
                        "prompt": "Make a diagnosis for 'infectious' or 'non-infectious'. If there are few taxa called in the below threshold data, evaluate each whether they are a candidate for contamination or classification artifact ('non-infectious'), or whether they are a candidate for a pathogen infectious'. Consider taxa even if they are rare human pathogens. If known viral human pathogens are called prioritize for 'infectious'",
                        "instructions": "Answer 'yes' if the integrated below threshold and target list data supports infectious diagnosis; 'no' if it does not or if no taxa were called."
                },
                "check": "llm_eval",
                "true_node": "diagnose_infectious",
                "false_node": "diagnose_non_infectious_flight_check"
            },
            "diagnose_infectious": {
                "question": {
                        "prompt": "Based on the synthesis of metagenomics evidence select the most likely pathogen candidates considering medical signficance (even of rare pathogens), metagenomics assay, sources of contamination, patient clinical history and sample type. Consider the contaminating taxa in your selection - report only the pathogen candidates with good supporting evidence for human infection in the clinical context provided..",
                        "instructions": "Select the most likely pathogen candidates from the data. DO NOT RETURN EXPLANATIONS, YOU MUST RETURN ONLY THE NAME OF THE CANDIDATE PATHOGENS IN GTDB FORMAT AS A STRING: 's__{Genus species}'."
                },
                "check": "llm_eval",
                "next": "describe_infectious"
            },
            "describe_infectious": {
                "question": {
                        "prompt": "You have selected one or more taxa as pathogen candidates and have the clinical information for this sample.",
                        "instructions": "Provide one paragraph with a precise summary of the role of each taxon (or group of related candidate pathogen taxa) in the context of the patient clinical data and the strength of the metagenomic profiling evidence, as well as a short descriptor of its microbiological and clinical background (or absence thereof) in this sample type. Keep sentences minimal with all sufficient information."
                },
                "check": "llm_eval",
                "next": "select_infectious"
            },
            "select_infectious": {
                "question": {
                        "prompt": "You have selected on or more taxa as pathogen candidates. Take into consideration the strength and weaknesses of the evidence from metagenomics assay results if multiple pathogen candidates were chosen. When encountering multiple bacterial or eukaryotic pathogen candidates often the causative agent is the one with the strongest evidence or abundance from metagenomic data especially the one with most assembled bases. If a virus is present that is a human pathogen and is a reasonable candidate given the clinical patient and sample type, prefer it over other candidates.",
                        "instructions": "Select a single pathogen from the pathogen candidates as most likely cause of infection. DO NOT RETURN EXPLANATIONS, YOU MUST RETURN ONLY THE NAME OF THE CANDIDATE PATHOGENS IN GTDB FORMAT AS A STRING: 's__{Genus species}'"
                },
                "check": "llm_eval",
                "next": "describe_select_infectious"
            },
            "describe_select_infectious": {
                "question": {
                        "prompt": "You have selected a single taxon as pathogen candidate and have the clinical information for this sample.",
                        "instructions": "Provide one paragraph with a precise summary of the role of this taxon in the context of the patient clinical data and the strength of the metagenomic profiling evidence, as well as a short descriptor of its microbiological and clinical background (or absence thereof) in this sample type. Explain why this taxon was chosen over other candidates as most likely pathogen in the sample."
                },
                "check": "llm_eval",
                "final_node": true
            },
            "diagnose_non_infectious_flight_check": {
                "question": {
                        "prompt": "Consider the below threshold and target evidence and consider whether changing to an 'infectious' call especially if the taxon is human infection associated but rare or unusual; alternatively make a call for 'non-infectious' when no taxa qualify as pathogen candidates or no taxa were called.",
                        "instructions": "Answer 'yes' if unusual or rare taxa were called that support an 'infectious' diagnosis; 'no' if it does not or if no taxa were called."
                },
                "check": "flight_check_non_infectious",
                "true_node": "diagnose_infectious_review",
                "false_node": "diagnose_non_infectious"
            },
            "diagnose_non_infectious": {
                "question": {
                        "prompt": "You have diagnosed this sample as non-infectious and have the clinical information for this sample. ",
                        "instructions": "Provide a short paragraph with a precise summary of the reasons for calling this sample non-infectious, highlighting conclusions around pathogen and contamination, and differentiating it with diagnostic results from the decision trees for differential tumor or infectious diagnosis."
                },
                "check": "llm_eval",
                "final_node": true
            },
            "diagnose_tumor": {
                "final_node": true
            },
            "diagnose_infectious_review": {
                "question": {
                    "prompt": "Based on the synthesis of metagenomics evidence select the most likely pathogen candidates considering medical signficance (even of rare pathogens), metagenomics assay, sources of contamination, patient clinical history and sample type (CSF or vitreous fluid). Make sure to differentiate between the most likely pathogen candidates with supporting evidence and contaminating taxa in your selection - report only the pathogen candidates with good supporting evidence for human infection in the clinical context provided.",
                    "instructions": "Select the most likely pathogen candidates from the data. DO NOT RETURN EXPLANATIONS, YOU MUST RETURN ONLY THE NAME OF THE CANDIDATE PATHOGENS IN GTDB FORMAT AS A STRING: 's__{Genus species}'."
                },
                "check": "llm_eval",
                "next": "describe_infectious_review"
            },
            "describe_infectious_review": {
                "question": {
                        "prompt": "You have selected on or more taxa as pathogen candidates for review and have the clinical information for this sample.",
                        "instructions": "Provide one paragraph with a precise summary of the role of each taxon (or group of related candidate pathogen taxa) in the context of the patient clinical data and the strength of the metagenomic profiling evidence, as well as a short descriptor of its microbiological and clinical background (or absence thereof) in this sample type."
                },
                "check": "llm_eval",
                "next": "select_infectious_review"
            },
            "select_infectious_review": {
                "question": {
                        "prompt": "You have selected one or more taxa as pathogen candidates with review. Take into consideration the strength and weaknesses of the evidence from metagenomics assay results if multiple pathogen candidates were chosen. When encountering multiple bacterial or eukaryotic pathogen candidates often the causative agent is the one with the strongest evidence or abundance from metagenomic data especially the one with most assembled bases. If a virus is present that is a human pathogen and is a reasonable candidate given the clinical patient and sample type, prefer it over other candidates. If a candidate is well known for being a human pathogen and is less frequently observed as contamination, prefer it over a candidate that can be a human pathogen but is more frequently observed as contamination.",
                        "instructions": "Select a single pathogen from the pathogen candidates as most likely cause of infection. DO NOT RETURN EXPLANATIONS, YOU MUST RETURN ONLY THE NAME OF THE CANDIDATE PATHOGENS IN GTDB FORMAT AS A STRING: 's__{Genus species}'"
                },
                "check": "llm_eval",
                "next": "describe_select_infectious_review"
            },
            "describe_select_infectious_review": {
                "question": {
                        "prompt": "You have selected a single taxon as pathogen candidate and have the clinical information for this sample.",
                        "instructions": "Provide one paragraph with a precise summary of the role of this taxon in the context of the patient clinical data and the strength of the metagenomic profiling evidence, as well as a short descriptor of its microbiological and clinical background (or absence thereof) in this sample type. Explain why this taxon was chosen over other candidates as most likely pathogen in the sample."
                },
                "check": "llm_eval",
                "final_node": true
            }
        }
        "#;

        let tree: DecisionTree = serde_json::from_str(decision_tree_json)?;
        
        let llm_client = async_openai::Client::new();
        let state = AgentState::new();

        // Build a simple knowledge graph from the decision tree.
        let mut knowledge_graph = Graph::<String, String>::new();
        let mut node_indices = HashMap::new();

        // Create a node in the graph for each decision tree node.
        for (key, _node) in &tree {
            let idx = knowledge_graph.add_node(key.clone());
            node_indices.insert(key.clone(), idx);
        }
        // Add edges for transitions.
        for (key, node) in &tree {
            if let Some(ref true_node) = node.true_node {
                if let (Some(&from), Some(&to)) =
                    (node_indices.get(key), node_indices.get(true_node))
                {
                    knowledge_graph.add_edge(from, to, "true".to_string());
                }
            }
            if let Some(ref false_node) = node.false_node {
                if let (Some(&from), Some(&to)) =
                    (node_indices.get(key), node_indices.get(false_node))
                {
                    knowledge_graph.add_edge(from, to, "false".to_string());
                }
            }
            if let Some(ref next) = node.next {
                if let (Some(&from), Some(&to)) =
                    (node_indices.get(key), node_indices.get(next))
                {
                    knowledge_graph.add_edge(from, to, "next".to_string());
                }
            }
        }

        Ok(DiagnosticAgent {
            tree,
            cerebro_client,
            llm_client,
            state,
            knowledge_graph,
            model,
            diagnostic_memory,
            contam_history
        })
    }
    fn model_supports_system_message(&self) -> bool {
        Vec::from(["gpt-4o".to_string(), "gpt-4o-mini".to_string()]).contains(&self.model)
    }
    // Evaluate an LLM prompt using async-openai.
    async fn evaluate_diagnostic_llm(&self, question: &Question, context: &str) -> Result<bool> {
        let prompt_text = question.to_prompt();

        let messages = if self.model_supports_system_message() {
            vec![
                ChatCompletionRequestMessage::System("You are a diagnostic assistant. Answer 'yes' or 'no' only.".into()),
                ChatCompletionRequestMessage::User(format!("{}\n\nCase Context:\n{}", prompt_text, context).into()),
            ]
        } else {
            vec![
                ChatCompletionRequestMessage::User(format!("{}\n\nCase Context:\n{}", prompt_text, context).into()),
            ]
        };

        let request = CreateChatCompletionRequestArgs::default()
            .model(&self.model)
            .messages(messages)
            .build()?;
        let response = self.llm_client.chat().create(request).await?;
        let reply = &response.choices[0].message.content;
        Ok(reply.as_ref().map(|s| s.to_lowercase().contains("yes")).unwrap_or(false))
    }

    // Evaluate an LLM prompt using async-openai.
    async fn evaluate_llm(&self, question: &Question, context: &str) -> Result<String> {
        let prompt_text = question.to_prompt();
        
        let messages = if self.model_supports_system_message() {
            vec![
                ChatCompletionRequestMessage::System("You are a diagnostic assistant.".into()),
                ChatCompletionRequestMessage::User(format!("{}\n\nCase Context:\n{}", prompt_text, context).into()),
            ]
        } else {
            vec![
                ChatCompletionRequestMessage::User(format!("{}\n\nCase Context:\n{}", prompt_text, context).into()),
            ]
        };

        let request = CreateChatCompletionRequestArgs::default()
            .model(&self.model)
            .messages(messages)
            .build()?;
        let response = self.llm_client.chat().create(request).await?;
        let reply = &response.choices[0].message.content;
        Ok(reply.as_ref().unwrap_or(&format!("Unknown response from LLM")).to_string())
    }

    /// Primes the LLM with the background explanation.
    pub async fn prime_llm(&self, background: &str) -> Result<()> {
        let messages = if self.model_supports_system_message() {
            vec![
                ChatCompletionRequestMessage::System(background.into()),
                ChatCompletionRequestMessage::User("Please acknowledge that you have received the background context.".into()),
            ]
        } else {
            vec![
                ChatCompletionRequestMessage::User("Please acknowledge that you have received the background context.".into()),
            ]
        };
        let request = CreateChatCompletionRequestArgs::default()
            .model(&self.model)
            .messages(messages)
            .build()?;
        let _ = self.llm_client.chat().create(request).await?;
        Ok(())
    }
    /// Gets the current state memory into a JSON string.
    fn _get_agent_state_memory(&self) -> String {
        serde_json::to_string_pretty(&self.state.memory).unwrap_or_else(|_| "{}".into())
    }
    /// Gets the current diagnostic history into a JSON string.
    fn get_history(&self) -> String {
        serde_json::to_string_pretty(&self.state.history).unwrap_or_else(|_| "{}".into())
    }
    /// Gets the current threshold considerations
    fn _get_threshold_considerations(&self) -> String {
        let atc = self.state.memory.get("data__AboveThresholdQuery__Considerations").map_or("[]", |x| x.deref());
        let btc = self.state.memory.get("data__BelowThresholdQuery__Considerations").map_or("[]", |x| x.deref());
        let ttc = self.state.memory.get("data__TargetThresholdQuery__Considerations").map_or("[]", |x| x.deref());
        format!("Above threshold filter configuration: {atc} === Below threshold filter configuration: {btc} === Target threshold filter contamination {ttc}")
    }
    /// Gets the current threshold candidate taxa by deserializing JSON strings to Vec<Taxon>
    fn get_threshold_candidates(&self) -> ThresholdCandidates {

        // Retrieve JSON strings from memory, defaulting to "[]" if not found.
        let atc_json = self.state.memory.get("data__AboveThresholdQuery__CandidateTaxa").map_or("[]", |x| x.deref());
        let btc_json = self.state.memory.get("data__BelowThresholdQuery__CandidateTaxa").map_or("[]", |x| x.deref());
        let ttc_json = self.state.memory.get("data__TargetThresholdQuery__CandidateTaxa").map_or("[]", |x| x.deref());
        
        ThresholdCandidates { 
            above_threshold: serde_json::from_str(atc_json).ok(), 
            below_threshold: serde_json::from_str(btc_json).ok(), 
            target_list: serde_json::from_str(ttc_json).ok() 
        }
    }
    fn threshold_candidates_to_str(&self, full_context: &mut String, threshold_candidates: &ThresholdCandidates) {
        full_context.push_str(&threshold_candidates.to_str());
    }
    /// Runs the diagnostic pathway and returns a DiagnosticResult as JSON.
    pub async fn run(&mut self, gp_config: &mut MetaGpConfig, clinical_context: &str) -> Result<DiagnosticResult> {

        // Create the identifier query for this sample:
        self.state.log("identifier_schema", &gp_config.identifiers.to_json());

        let mut node_key = "start".to_string();

        // Loop through the decision tree.
        while let Some(node) = self.tree.get(&node_key) {
            log::info!("Processing node: {}", node_key);

            self.state.log("current_node", &node_key);
            
            // Record the current node in history (include question if available)
            let question_text = node.question.as_ref().map(|q| q.to_prompt());

            match node.check {
                Some(CheckType::LlmEval) | Some(CheckType::LlmDiagnosticEval) | Some(CheckType::LlmJsonEval) => {
                    
                    let question = node.question.as_ref().expect("LLM JSON eval requires a question");
                    
                    // Special handling for the output_evaluation node:
                    if node_key == "output_evaluation" {

                        // Get JSON-structured memory injection.
                        let diagnostic_memory = self.get_history();

                        let threshold_candidate_taxa = self.get_threshold_candidates();

                        // Combine clinical context with the injected memory.
                        let mut full_context = format!("Clinical Context: {}\n\n{}", clinical_context,  if self.diagnostic_memory { format!("Diagnostic memory: {}\n\n", diagnostic_memory) } else { String::new() });

                        // Combine data context with the query
                        self.threshold_candidates_to_str(&mut full_context, &threshold_candidate_taxa);

                        // Fetch JSON output from the LLM with the full context.
                        let diagnostic_response = self.evaluate_diagnostic_llm(question, &full_context).await?;
                        
                        // Log the JSON response in history.
                        self.state.add_history(&node_key, Some(question.to_prompt()), Some(&diagnostic_response.to_string()));

                        if diagnostic_response {
                            node_key = node.true_node.clone().expect("True node key expected");
                        } else {
                            node_key = node.false_node.clone().expect("False node key expected");
                        }

                    } else if node_key == "diagnose_infectious" {

                        // Get JSON-structured memory injection.
                        let diagnostic_memory = self.get_history();

                        let threshold_candidate_taxa = self.get_threshold_candidates();

                        // Combine clinical context with the injected memory.
                        let mut full_context = format!("Clinical Context: {}\n\n{}", clinical_context, if self.diagnostic_memory { format!("Diagnostic memory: {}\n\n", diagnostic_memory) } else { String::new() });

                        // Combine data context with the query
                        self.threshold_candidates_to_str(&mut full_context, &threshold_candidate_taxa);

                        // Fetch JSON output from the LLM with the full context.
                        let pathogen_candidates = self.evaluate_llm(question, &full_context).await?;
                        
                        // Log the JSON response in history.
                        self.state.add_history(&node_key, Some(question.to_prompt()), Some(&pathogen_candidates.to_string()));

                        // Save the parsed evaluation details in memory.
                        self.state.log("pathogen_candidates", &pathogen_candidates.to_string());
                        
                        if let Some(next_node) = &node.next {
                            node_key = next_node.to_string();
                        }

                        
                    } else if node_key == "describe_infectious" {

                        // Get JSON-structured memory injection.
                        let diagnostic_memory = self.get_history();

                        let pathogen_candidates = self.state.memory.get("pathogen_candidates");

                        // Combine clinical context with the injected memory.
                        let full_context = format!("Clinical Context: {}\n\n{}\n\nPathogen candidate taxa: {:?}", clinical_context, if self.diagnostic_memory { format!("Diagnostic memory: {}", diagnostic_memory) } else { String::new() }, pathogen_candidates);

                        // Fetch JSON output from the LLM with the full context.
                        let pathogen_description = self.evaluate_llm(question, &full_context).await?;
                        
                        // Log the JSON response in history.
                        self.state.add_history(&node_key, Some(question.to_prompt()), Some(&pathogen_description.to_string()));

                        // Save the parsed evaluation details in memory.
                        self.state.log("pathogen_candidate_description", &pathogen_description.to_string());

                        if let Some(next_node) = &node.next {
                            node_key = next_node.to_string();
                        }
                        
                    } else if node_key == "select_infectious" {

                        // Get JSON-structured memory injection.
                        let diagnostic_memory = self.get_history();

                        let threshhold_candidates = self.get_threshold_candidates();
                        let pathogen_candidates = self.state.memory.get("pathogen_candidates");
                        let pathogen_candidate_description = self.state.memory.get("pathogen_candidate_description");

                        // Combine clinical context with the injected memory.
                        let full_context = format!("Clinical Context: {}\n\n{}\n\nPathogen candidate taxa: {:?}\n\nPathogen candidate description: {:?}\n\nThreshold candidates: {}", clinical_context, if self.diagnostic_memory { format!("Diagnostic memory: {}", diagnostic_memory) } else { String::new() }, pathogen_candidates, pathogen_candidate_description, threshhold_candidates.to_str());

                        // Fetch JSON output from the LLM with the full context.
                        let pathogen_selection = self.evaluate_llm(question, &full_context).await?;
                        
                        // Log the JSON response in history.
                        self.state.add_history(&node_key, Some(question.to_prompt()), Some(&pathogen_selection));

                        // Save the parsed evaluation details in memory.
                        self.state.log("pathogen_selection", &pathogen_selection);

                        if let Some(next_node) = &node.next {
                            node_key = next_node.to_string();
                        }
                        
                    } else if node_key == "describe_select_infectious" {

                        // Get JSON-structured memory injection.
                        let diagnostic_memory = self.get_history();

                        let pathogen_selection = self.state.memory.get("pathogen_selection");

                        // Combine clinical context with the injected memory.
                        let full_context = format!("Clinical Context: {}\n\n{}\n\nSelected pathogen candidate: {:?}", clinical_context, if self.diagnostic_memory { format!("Diagnostic memory: {}", diagnostic_memory) } else { String::new() }, pathogen_selection);

                        // Fetch JSON output from the LLM with the full context.
                        let pathogen_description = self.evaluate_llm(question, &full_context).await?;
                        
                        // Log the JSON response in history.
                        self.state.add_history(&node_key, Some(question.to_prompt()), Some(&pathogen_description));

                        // Save the parsed evaluation details in memory.
                        self.state.log("pathogen_description", &pathogen_description);
                        
                    } else if node_key == "integrate_below_target_evidence" {

                        // Get JSON-structured memory injection.
                        let diagnostic_memory = self.get_history();

                        let threshold_candidate_taxa = self.get_threshold_candidates();

                        // Combine clinical context with the injected memory.
                        let mut full_context = format!("Clinical Context: {}\n\n{}\n\n", clinical_context, if self.diagnostic_memory { format!("Diagnostic memory: {}", diagnostic_memory) } else { String::new() });

                        // Combine data context with the query
                        self.threshold_candidates_to_str(&mut full_context, &threshold_candidate_taxa);

                        // Fetch JSON output from the LLM with the full context.
                        let diagnostic_response = self.evaluate_diagnostic_llm(question, &full_context).await?;
                        
                        // Log the JSON response in history.
                        self.state.add_history(&node_key, Some(question.to_prompt()), Some(&diagnostic_response.to_string()));

                        if diagnostic_response {
                            node_key = node.true_node.clone().expect("True node key expected");
                        } else {
                            node_key = node.false_node.clone().expect("False node key expected");
                        }
                        
                    } else if node_key == "diagnose_infectious_review" {

                        // Get JSON-structured memory injection.
                        let diagnostic_memory = self.get_history();

                        let threshold_candidate_taxa = self.get_threshold_candidates();

                        // Combine clinical context with the injected memory.
                        let mut full_context = format!("Clinical Context: {}\n\n{}\n\n", clinical_context, if self.diagnostic_memory { format!("Diagnostic memory: {}", diagnostic_memory) } else { String::new() });
                        
                        // Combine data context with the query
                        self.threshold_candidates_to_str(&mut full_context, &threshold_candidate_taxa);

                        // Fetch JSON output from the LLM with the full context.
                        let pathogen_candidates = self.evaluate_llm(question, &full_context).await?;
                        
                        // Log the JSON response in history.
                        self.state.add_history(&node_key, Some(question.to_prompt()), Some(&pathogen_candidates.to_string()));

                        // Save the parsed evaluation details in memory.
                        self.state.log("pathogen_candidates", &pathogen_candidates.to_string());
                        
                        if let Some(next_node) = &node.next {
                            node_key = next_node.to_string();
                        }
                        
                    } else if node_key == "describe_infectious_review" {

                        // Get JSON-structured memory injection.
                        let diagnostic_memory = self.get_history();

                        let pathogen_candidates = self.state.memory.get("pathogen_candidates");

                        // Combine clinical context with the injected memory.
                        let full_context = format!("Clinical Context: {}\n\n{}\n\nPathogen candidate taxa: {:?}", clinical_context, if self.diagnostic_memory { format!("Diagnostic memory: {}", diagnostic_memory) } else { String::new() }, pathogen_candidates);

                        // Fetch JSON output from the LLM with the full context.
                        let pathogen_description = self.evaluate_llm(question, &full_context).await?;
                        
                        // Log the JSON response in history.
                        self.state.add_history(&node_key, Some(question.to_prompt()), Some(&pathogen_description.to_string()));

                        // Save the parsed evaluation details in memory.
                        self.state.log("pathogen_candidate_description", &pathogen_description.to_string());

                        if let Some(next_node) = &node.next {
                            node_key = next_node.to_string();
                        }
                        
                    } else if node_key == "select_infectious_review" {

                        // Get JSON-structured memory injection.
                        let diagnostic_memory = self.get_history();

                        let threshhold_candidates = self.get_threshold_candidates();
                        let pathogen_candidates = self.state.memory.get("pathogen_candidates");
                        let pathogen_candidate_description = self.state.memory.get("pathogen_candidate_description");

                        // Combine clinical context with the injected memory.
                        let full_context = format!("Clinical Context: {}\n\n{}\n\nPathogen candidate taxa: {:?}\n\nPathogen candidate description: {:?}\n\nThreshold candidate taxa: {}", clinical_context, if self.diagnostic_memory { format!("Diagnostic memory: {}", diagnostic_memory) } else { String::new() }, pathogen_candidates, pathogen_candidate_description, threshhold_candidates.to_str());

                        // Fetch JSON output from the LLM with the full context.
                        let pathogen_selection = self.evaluate_llm(question, &full_context).await?;
                        
                        // Log the JSON response in history.
                        self.state.add_history(&node_key, Some(question.to_prompt()), Some(&pathogen_selection.to_string()));

                        // Save the parsed evaluation details in memory.
                        self.state.log("pathogen_selection", &pathogen_selection.to_string());

                        if let Some(next_node) = &node.next {
                            node_key = next_node.to_string();
                        }
                        
                    } else if node_key == "describe_select_infectious_review" {

                        // Get JSON-structured memory injection.
                        let diagnostic_memory = self.get_history();

                        let pathogen_selection = self.state.memory.get("pathogen_selection");

                        // Combine clinical context with the injected memory.
                        let full_context = format!("Clinical Context: {}\n\n{}\n\nSelected pathogen candidate: {:?}", clinical_context, if self.diagnostic_memory { format!("Diagnostic memory: {}", diagnostic_memory) } else { String::new() }, pathogen_selection);

                        // Fetch JSON output from the LLM with the full context.
                        let pathogen_description = self.evaluate_llm(question, &full_context).await?;
                        
                        // Log the JSON response in history.
                        self.state.add_history(&node_key, Some(question.to_prompt()), Some(&pathogen_description));

                        // Save the parsed evaluation details in memory.
                        self.state.log("pathogen_description", &pathogen_description);
                        
                    } else {
                        // Otherwise full diagnostic decision LLM evaluation returning yes/no.
                        let answer = self.evaluate_llm(question, clinical_context).await?;

                        self.state.log("llm_decision", &answer);
                        self.state.add_history(&node_key, Some(question.to_prompt()), Some(&answer));
                        
                        if let Some(next_node) = &node.next {
                            node_key = next_node.to_string();
                        } else {
                            if node.final_node.unwrap_or(false) {
                                self.state.add_history(&node_key, question_text, Some("Final node reached"));
                                log::info!("Final node reached: {}.", node_key);
                                break;
                            }

                            if let Some(next_node) = &node.next {
                                node_key = next_node.to_string();
                            }
                        }
                    }
                },
                Some(CheckType::FlightCheckNonInfectious) => {

                    let question = node.question.as_ref().expect("LLM eval requires a question");

                    // Candidate taxa
                    let threshold_candidate_taxa = self.get_threshold_candidates();

                    // Get JSON-structured memory injection.
                    let above_threshold_answer = self.state.get_history(&"check_above_threshold").and_then(|d| d.answer.clone()).unwrap_or(String::from("no"));
                    let below_threshold_answer = self.state.get_history(&"check_below_threshold").and_then(|d| d.answer.clone()).unwrap_or(String::from("no"));
                    let target_threshold_answer = self.state.get_history(&"check_target_threshold").and_then(|d| d.answer.clone()).unwrap_or(String::from("no"));

                    // Combine clinical context with the injected memory.
                    let mut full_context = format!(
                        "Clinical Context: {}\n\nAbove threshold answer for infectious sample: {:?}\n\nBelow threshold answer for infectious sample: {:?}\n\nTarget list answer for infectious sample: {:?}", 
                        clinical_context, above_threshold_answer, below_threshold_answer, target_threshold_answer
                    );
                    
                    // Combine data context with the query
                    self.threshold_candidates_to_str(&mut full_context, &threshold_candidate_taxa);

                    // Fetch JSON output from the LLM with the full context.
                    let diagnostic_result = self.evaluate_diagnostic_llm(question, &full_context).await?;
                    
                    // Log the JSON response in history.
                    self.state.add_history(&node_key, Some(question.to_prompt()), Some(&diagnostic_result.to_string()));

                    // Save the parsed evaluation details in memory.
                    self.state.log("non_infectious_check", &diagnostic_result.to_string());
                    
                    if let Some(next_node) = &node.next {
                        node_key = next_node.to_string();
                    } else {
                        node_key = if diagnostic_result {
                            node.true_node.clone()
                        } else {
                            node.false_node.clone()
                        }
                        .expect("Next node key expected");
                    }

                }
                Some(CheckType::AboveThresholdQuery) => {

                    let (above_threshold_taxa, above_threshold_contam) = self.cerebro_client.get_taxa( 
                        &CerebroIdentifierSchema::from_gp_config(gp_config), 
                        &TaxonFilterConfig::gp_above_threshold(gp_config.ignore_taxstr.clone()), 
                        &mut gp_config.contamination,
                        self.contam_history
                    )?;

                    self.state.log(
                        "data__AboveThresholdQuery__CandidateTaxa", 
                        &serde_json::to_string_pretty(&above_threshold_taxa).unwrap_or("{}".to_string())
                    );
                    self.state.log(
                        "data__AboveThresholdQuery__ContamTaxa", 
                        &serde_json::to_string_pretty(&above_threshold_contam).unwrap_or("{}".to_string())
                    );
                    
                    // Combine data context with the query
                    let mut full_context = String::new();
                    
                    // Combine data context with the query
                    self.threshold_candidates_to_str(&mut full_context, &ThresholdCandidates::from_above_threshold(above_threshold_taxa));

                    let question = node.question.as_ref().expect("LLM diagnostic eval requires a question");

                    // Fetch bool output from the LLM with the full context.
                    let diagnostic_response = self.evaluate_diagnostic_llm(question, &full_context).await?;

                    self.state.add_history(&node_key, node.question.as_ref().map(|q| q.to_prompt()), Some(&diagnostic_response.to_string()));

                    if let Some(next_node) = &node.next {
                        node_key = next_node.to_string();
                    } else {
                        node_key = if diagnostic_response {
                            node.true_node.clone()
                        } else {
                            node.false_node.clone()
                        }
                        .expect("Next node key expected");
                    }

                }
                Some(CheckType::BelowThresholdQuery) => {
                    
                    let (below_threshold_taxa, below_threshold_contam) = self.cerebro_client.get_taxa(
                        &CerebroIdentifierSchema::from_gp_config(gp_config), 
                        &TaxonFilterConfig::gp_below_threshold(gp_config.ignore_taxstr.clone()), 
                        &mut gp_config.contamination,
                        false
                    )?;

                    self.state.log(
                        "data__BelowThresholdQuery__CandidateTaxa", 
                        &serde_json::to_string_pretty(&below_threshold_taxa).unwrap_or("{}".to_string())
                    );
                    self.state.log(
                        "data__BelowThresholdQuery__ContamTaxa", 
                        &serde_json::to_string_pretty(&below_threshold_contam).unwrap_or("{}".to_string())
                    );


                    // Combine data context with the query
                    let mut full_context = String::new();
                    
                    // Combine data context with the query
                    self.threshold_candidates_to_str(&mut full_context, &ThresholdCandidates::from_below_threshold(below_threshold_taxa));

                    let question = node.question.as_ref().expect("LLM diagnostic eval requires a question");

                    // Fetch bool output from the LLM with the full context.
                    let diagnostic_response = self.evaluate_diagnostic_llm(question, &full_context).await?;

                    self.state.add_history(&node_key, node.question.as_ref().map(|q| q.to_prompt()), Some(&diagnostic_response.to_string()));

                    if let Some(next_node) = &node.next {
                        node_key = next_node.to_string();
                    } else {
                        node_key = if diagnostic_response {
                            node.true_node.clone()
                        } else {
                            node.false_node.clone()
                        }
                        .expect("Next node key expected");
                    }
                }
                Some(CheckType::TargetThresholdQuery) => {

                    let (target_threshold_taxa, target_threshold_contam) = self.cerebro_client.get_taxa(
                        &CerebroIdentifierSchema::from_gp_config(gp_config), 
                        &TaxonFilterConfig::gp_target_threshold(gp_config.ignore_taxstr.clone()), 
                        &mut gp_config.contamination,
                        false
                    )?;
            
                    self.state.log(
                        "data__TargetThresholdQuery__CandidateTaxa", 
                        &serde_json::to_string_pretty(&target_threshold_taxa).unwrap_or("{}".to_string())
                    );
                    self.state.log(
                        "data__TargetThresholdQuery__ContamTaxa", 
                        &serde_json::to_string_pretty(&target_threshold_contam).unwrap_or("{}".to_string())
                    );

                    // Combine data context with the query
                    let mut full_context = String::new();

                    // Combine data context with the query
                    self.threshold_candidates_to_str(&mut full_context, &ThresholdCandidates::from_target_list(target_threshold_taxa));

                    let question = node.question.as_ref().expect("LLM diagnostic eval requires a question");

                    // Fetch bool output from the LLM with the full context.
                    let diagnostic_response = self.evaluate_diagnostic_llm(question, &full_context).await?;

                    self.state.add_history(&node_key, node.question.as_ref().map(|q| q.to_prompt()), Some(&diagnostic_response.to_string()));

                    if let Some(next_node) = &node.next {
                        node_key = next_node.to_string();
                    } else {
                        node_key = if diagnostic_response {
                            node.true_node.clone()
                        } else {
                            node.false_node.clone()
                        }
                        .expect("Next node key expected");
                    }

                },
                Some(CheckType::AneuploidyQuery) => {

                    let aneuploidy_data = self.cerebro_client.get_aneuploidy(&gp_config.sample)?;
                    let considerations = "[Host aneuploidy evaluation from host genome-wide segmental copy number variation (CNV) - large segmental aberrations across multiple chromosomes with 
                    data not widely spread around the genome-wide CNV estimate. If data is widely spread and diffuse it indicates insufficient host genome coverage. Positive large 
                    segmental abberrations in CNV across multiple chromosomes indicates tumorous host DNA.]";
                    
                    self.state.log(
                        "data__AneuploidyQuery__Considerations",
                        considerations 
                    );

                    self.state.log(
                        "data__AneuploidyQuery__Data", 
                        aneuploidy_data.unwrap_or("No data available")
                    );
                    
                    // Combine data context with the query
                    let full_context = format!("Aneuploidy data: {:#?}\n\nConsiderations: {}", aneuploidy_data, considerations);

                    let question = node.question.as_ref().expect("LLM diagnostic eval requires a question");

                    // Fetch bool output from the LLM with the full context.
                    let diagnostic_response = self.evaluate_diagnostic_llm(question, &full_context).await?;

                    self.state.add_history(&node_key, node.question.as_ref().map(|q| q.to_prompt()), Some(&diagnostic_response.to_string()));
                    
                    if let Some(next_node) = &node.next {
                        node_key = next_node.to_string();
                    } else {
                        node_key = if diagnostic_response {
                            node.true_node.clone()
                        } else {
                            node.false_node.clone()
                        }
                        .expect("Next node key expected");
                    }
                },
                Some(CheckType::TaxaHistoryQuery) => {

                    // let above_threshold_contam_taxa: Vec<Taxon> = serde_json::from_str(self.state.memory.get("data__AboveThresholdQuery__ContamTaxa").map_or("[]", |r| r.deref()))?;

                    // let taxa_history_data = self.cerebro_client.get_taxa_history(above_threshold_contam_taxa)?;

                    // let considerations = "[Historical occurrence of taxa in above threshold calls (sum of RPM from k-mer and alignment profiling) for DNA libraries as a function of human host background RPM - 
                    // positive linear regression of log-transformed data indicates problems with human sequence contamination in reference databases, negative values indicate normal expected 
                    // behaviour with increasing host biomass background to have reduced taxon calls; outliers from the regression 95% CI interval indicate real taxon calls and may negate the 
                    // identification of this taxon as a contaminant based on the prevalence contamination filter and occurrence across the majority of samples assessed in that filter]";

                    // self.state.log(
                    //     "data__TaxaHistoryData__Considerations", 
                    //     considerations
                    // );

                    // self.state.log(
                    //     "data__TaxaHistoryData__Data",
                    //     taxa_history_data.unwrap_or("No data available") 
                    // );

                    // // Combine data context with the query
                    // let full_context = format!("Taxa history data: {:#?}\n\nConsiderations: {}", taxa_history_data, considerations);

                    // let question = node.question.as_ref().expect("LLM diagnostic eval requires a question");

                    // // Fetch bool output from the LLM with the full context.
                    // let diagnostic_response = self.evaluate_diagnostic_llm(question, &full_context).await?;

                    // self.state.add_history(&node_key, node.question.as_ref().map(|q| q.to_prompt()), Some(&diagnostic_response.to_string()));
                    
                    
                    // if let Some(next_node) = &node.next {
                    //     node_key = next_node.to_string();
                    // } else {
                    //     node_key = if diagnostic_response {
                    //         node.true_node.clone()
                    //     } else {
                    //         node.false_node.clone()
                    //     }
                    //     .expect("Next node key expected");
                    // }
                }
                None => {},
            }

            // If this is a final node, log and exit the loop.
            if node.final_node.unwrap_or(false) {
                self.state.add_history(&node_key, question_text, Some("Final node reached"));
                log::info!("Final node reached: {}.", node_key);
                break;
            }
        };


        // Map the final node key to an overall diagnosis
        let diagnosis = match node_key.as_str() {
            "diagnose_infectious" => Diagnosis::Infectious,
            "describe_select_infectious" => Diagnosis::Infectious,
            "diagnose_infectious_review" => Diagnosis::InfectiousReview,
            "describe_select_infectious_review" => Diagnosis::InfectiousReview,
            "diagnose_tumor" => Diagnosis::Tumor,
            "diagnose_non_infectious" => Diagnosis::NonInfectious,
            "diagnose_non_infectious_review" => Diagnosis::NonInfectiousReview,
            _ => Diagnosis::Unknown,
        };

        let candidate = self.state.memory.get("pathogen_candidates").cloned();
        let pathogen = self.state.memory.get("pathogen_selection").cloned();

        let reason_non_infectious = self.state.memory.get("non_infectious_reason").cloned();
        let candidate_description = self.state.memory.get("pathogen_candidate_description").cloned();
        let pathogen_description = self.state.memory.get("pathogen_description").cloned();

        let taxa = ThresholdTaxa::from_agent_memory(&self.state.memory)?;

        let diagnostic_result = DiagnosticResult {
            diagnosis,
            candidate,
            candidate_description,
            pathogen,
            pathogen_description,
            reason_non_infectious,
            clinical_context: clinical_context.to_string(),
            diagnostic_history_shorthand: self.state.shorthand.clone(),
            diagnostic_history_annotated: self.state.history.clone(),
            taxa
        };
        

        log::info!(
            "Diagnosis: '{}'",
            format!("{}", diagnostic_result.diagnosis).blue()
        );

        if diagnostic_result.diagnosis == Diagnosis::Infectious || diagnostic_result.diagnosis == Diagnosis::InfectiousReview {
            log::info!(
                "Candidates: {} | Description: {}", 
                diagnostic_result.candidate.clone().unwrap_or_default().red(), 
                diagnostic_result.candidate_description.clone().unwrap_or_default()
            );        
            log::info!(
                "Pathogen: {} | Description: {}", 
                diagnostic_result.pathogen.clone().unwrap_or_default().red(), 
                diagnostic_result.pathogen_description.clone().unwrap_or_default()
            )
        }

        if diagnostic_result.diagnosis == Diagnosis::NonInfectious || diagnostic_result.diagnosis == Diagnosis::NonInfectiousReview {
            log::info!(
                "Reason: {}", 
                diagnostic_result.reason_non_infectious.clone().unwrap_or_default() 
            );     
        }

        Ok(diagnostic_result)
    }

    pub fn print_knowledge_graph(&self) {
        println!("Knowledge Graph:");
        for node in self.knowledge_graph.node_indices() {
            let label = &self.knowledge_graph[node];
            println!("Node {}: {}", node.index(), label);
        }
        for edge in self.knowledge_graph.edge_references() {
            let source = &self.knowledge_graph[edge.source()];
            let target = &self.knowledge_graph[edge.target()];
            println!("Edge: {} --({})--> {}", source, edge.weight(), target);
        }
    }

    pub fn print_decision_tree(&self, tree: &DecisionTree, root: &str) {

        // Recursive helper function
        fn recurse(tree: &DecisionTree, node_key: &str, prefix: &str, is_last: bool) {
            let branch = if is_last { "└── " } else { "├── " };
            // Print the current node key in bold white.
            println!("{}{}{}", prefix, branch, node_key.bold().white());
            
            // Prepare the prefix for the children.
            let new_prefix = if is_last { format!("{}    ", prefix) } else { format!("{}│   ", prefix) };

            // Collect children from the node: we assume true_node, false_node, and next.
            let mut children: Vec<(&str, &String)> = Vec::new();
            if let Some(node) = tree.get(node_key) {
                if let Some(child) = &node.true_node {
                    children.push(("true", child));
                }
                if let Some(child) = &node.false_node {
                    children.push(("false", child));
                }
                if let Some(child) = &node.next {
                    children.push(("next", child));
                }
            }

            let count = children.len();
            for (i, (label, child_key)) in children.into_iter().enumerate() {
                
                // For extra style, we can color the label (e.g., yellow).
                let colored_label = format!("{}: ", label).yellow();
                let colored_child = child_key.bold().white();
                let line = format!("{}{}", colored_label, colored_child);

                // Instead of printing the label separately, we can include it in the recursion:
                // Print the child node line.
                let is_last_child = i == count - 1;
                println!("{}{}{}", new_prefix, if is_last_child { "└── " } else { "├── " }, line);

                // Then recurse on this child if it exists in the tree.
                if tree.contains_key(child_key) {
                    // Adjust prefix for further children.
                    let next_prefix = if is_last_child { format!("{}    ", new_prefix) } else { format!("{}│   ", new_prefix) };
                    recurse(tree, child_key, &next_prefix, true);
                }
            }
        }

        // Print header
        println!("{}", "Decision Tree:".bold().underline());

        // Start recursion from the root.
        recurse(tree, root, "", true);
    }
}
