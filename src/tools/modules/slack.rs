use std::process;

use crate::tools::error::ToolError;


use actix_web_httpauth::headers::authorization::Bearer;
use reqwest::header::AUTHORIZATION;
use serde::{Serialize, Deserialize};


#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SlackMessage {
    pub channel: String,
    pub icon_emoji: String,
    pub blocks: Vec<SlackMessageBlock>
}

impl SlackMessage {
    pub fn new(channel: &str, message: &str) -> Self {
        Self {
            channel: channel.to_string(),
            icon_emoji: String::from(":satellite:"),
            blocks: vec![
            SlackMessageBlock::rich_text_block(
                vec![RichTextElement::from(
                    vec![RichTextObject::from(message)]
                )]
            )]
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SlackMessageBlock {
    pub r#type: String,
    pub elements: Vec<RichTextElement>,
    pub block_id: String,
}
impl SlackMessageBlock {
    pub fn rich_text_block(elements: Vec<RichTextElement>) -> Self {
        Self {
            r#type: String::from("rich_text"), elements, block_id: uuid::Uuid::new_v4().to_string()
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RichTextElement {
    pub r#type: String,
    pub elements: Vec<RichTextObject>,
}
impl RichTextElement {
    pub fn from(elements: Vec<RichTextObject>) -> Self {
        Self {
            r#type: String::from("rich_text_section"), elements
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RichTextObject {
    pub r#type: String,
    pub text: String
}
impl RichTextObject {
    pub fn from(text: &str) -> Self {
        Self {
            r#type: String::from("text"), text: text.to_string()
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SlackMessengerResponse {
    pub ok: bool
}

#[derive(Debug, Clone)]
pub struct SlackMessenger {
    pub url: String,
    pub token: String,
    pub client: reqwest::blocking::Client
}
impl SlackMessenger {
    pub fn new(token: &str) -> Self {
        let client = reqwest::blocking::Client::builder().build().expect("Failed to initialize Reqwest blocking client!");
        Self { url: String::from("https://slack.com/api/chat.postMessage"), token: token.to_string(), client }
    }
    pub fn send(&self, message: &SlackMessage) -> Result<(), ToolError> {

        let response: reqwest::blocking::Response = self.client.post(self.url.clone())
            .header(AUTHORIZATION, Bearer::new(self.token.clone()).to_string())
            .json(&message)
            .send().expect("Failed to send POST request to Slack API");

        let response_data: SlackMessengerResponse = response.json().expect("Failed to read JSON response from Slack API");

        if !response_data.ok {
            log::error!("Failed to send notification to channel: {}", message.channel);
            process::exit(1)
        } else {
            log::info!("Successfully sent notification to channel: {}", message.channel);
        }

        Ok(())

    }
}
