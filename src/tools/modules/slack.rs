use std::process;
use crate::tools::error::ToolError;
use actix_web_httpauth::headers::authorization::Bearer;
use reqwest::header::AUTHORIZATION;
use serde::{Serialize, Deserialize};



#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SlackMessage {
    pub channel: String,
    pub icon_emoji: String,
    pub blocks: Vec<SlackMessageSectionBlock>
}

impl SlackMessage {
    pub fn new(channel: &str, message: &str) -> Self {
        Self {
            channel: channel.to_string(),
            icon_emoji: String::from(":satellite:"),
            blocks: vec![
                SlackMessageSectionBlock::new(TextObject::markdown(message))
            ]
        }
    }
    pub fn from(channel: &str, blocks: Vec<SlackMessageSectionBlock>) -> Self {
        Self {
            channel: channel.to_string(),
            icon_emoji: String::from(":satellite:"),
            blocks
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TextObject {
    pub r#type: String,
    pub text: String
}
impl TextObject {
    pub fn markdown(text: &str) -> Self {
        Self {
            r#type: String::from("mrkdwn"), text: text.to_string()
        }
    }
}


#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SlackMessageSectionBlock {
    pub r#type: String,
    pub text: TextObject,
    pub block_id: String,
}
impl SlackMessageSectionBlock {
    pub fn new(text: TextObject) -> Self {
        Self {
            r#type: String::from("section"), text, block_id: uuid::Uuid::new_v4().to_string()
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
        }

        Ok(())

    }
}
