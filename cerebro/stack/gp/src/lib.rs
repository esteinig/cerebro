pub mod error;
pub mod terminal;
pub mod utils;
pub mod gpt;

#[cfg(feature = "local")]
pub mod llama;
#[cfg(feature = "local")]
pub mod quantized;
#[cfg(feature = "local")]
pub mod qwen;
#[cfg(feature = "local")]
pub mod text;