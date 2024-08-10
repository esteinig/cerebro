use std::process;
use actix_web_httpauth::headers::authorization::Bearer;
use reqwest::header::AUTHORIZATION;
use serde::{Serialize, Deserialize};
use crate::{error::WatcherError, terminal::WatchArgs};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SlackConfig {
    pub channel: String,
    pub token: String
}
impl SlackConfig {
    pub fn from_args(watch_args: &WatchArgs) -> Option<Self> {
        match (&watch_args.slack_channel, &watch_args.slack_token) {
            (Some(channel), Some(token)) => Some(
                SlackConfig { channel: channel.to_string(), token: token.to_string() }
            ),
            _ => {
                log::info!("Slack channel and token configuration not provided for this watcher");
                None
            }
        }
    }
}

#[derive(Debug, Clone)]
pub struct SlackTools {
    pub client: SlackClient,
    pub message: MessageGenerator
}
impl SlackTools {
    pub fn from_config(slack_config: &SlackConfig) -> Self {
        Self {
            client: SlackClient::new(&slack_config.token),
            message: MessageGenerator::new(&slack_config.channel)
        }
    }
}

#[derive(Debug, Clone)]
pub struct MessageGenerator {
    pub channel: String
}
impl MessageGenerator {
    pub fn new(channel: &str) -> Self {
        Self { channel: channel.to_string() }
    }
    pub fn watcher_setup(&self, watcher_name: &str, watcher_loc: &str) -> SlackMessage {
        SlackMessage::new(&self.channel, &format!("[{watcher_name}@{watcher_loc}] New watcher has been initialised"))
    }
    pub fn input_detected(&self, run_id: &str, watcher_name: &str, watcher_loc: &str) -> SlackMessage {
        SlackMessage::new(&self.channel, &format!("[{watcher_name}@{watcher_loc}::{run_id}] New run input detected"))
    }
    // pub fn input_validation(&self, validation: &InputValidation, run_id: &str, watcher_name: &str, watcher_loc: &str) -> SlackMessage {
    //     match validation.pass() {
    //         true => SlackMessage::new(
    //             &self.channel, 
    //             &format!("[{watcher_name}@{watcher_loc}::{run_id}] Input validation passed")
    //         ),
    //         false => {
    //             let msg = &format!("[{watcher_name}@{watcher_loc}::{run_id}] *Input validation failed*");
    //             SlackMessage::from(&self.channel, vec![
    //                 vec![SlackMessageSectionBlock::new(TextObject::markdown(msg))],
    //                 validation.get_message_blocks(&run_id)
    //             ].concat())
    //         }
    //     }
    // }
}

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
pub struct SlackClient {
    pub url: String,
    pub token: String,
    pub client: reqwest::blocking::Client
}
impl SlackClient {
    pub fn new(token: &str) -> Self {
        let client = reqwest::blocking::Client::builder().build().expect("Failed to initialize Reqwest blocking client!");
        Self { url: String::from("https://slack.com/api/chat.postMessage"), token: token.to_string(), client }
    }
    pub fn send(&self, message: &SlackMessage) -> Result<(), WatcherError> {

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
