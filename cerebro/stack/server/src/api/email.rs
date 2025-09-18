
use handlebars::Handlebars;
use lettre::{
    message::header::ContentType,
    transport::smtp::authentication::Credentials,
    AsyncSmtpTransport,
    AsyncTransport,
    Message,
    Tokio1Executor,
};

use cerebro_model::api::config::Config;
use cerebro_model::api::users::model::User;

pub struct TemplateData {
    subject: String,
    name: String,
    msg: String,
    url: String,
    gpg: String,
}
impl TemplateData {
    pub fn new(subject: &str, name: &str, msg: &str, url: &str, gpg: &str) -> Self {
        Self { subject: subject.to_owned(), name: name.to_owned(), msg: msg.to_owned(), url: url.to_owned(), gpg: gpg.to_owned()}
    }
    pub fn to_json(&self) -> serde_json::Value {
        serde_json::json!({
            "subject": &self.subject,
            "name": &self.name,
            "msg": &self.msg,
            "url": &self.url,
            "gpg": &self.gpg
        })
    }
}

pub struct Email {
    pub user: User,
    pub from: String,
    pub config: Config,
}
impl Email {
    pub fn new(user: &User, config: &Config) -> Self {
        Email { 
            user: user.to_owned(),
            from: format!("{} <{}>", config.smtp.from.to_owned(), config.smtp.username.to_owned()), 
            config: config.to_owned(),
        }
    }

    fn new_transport(
        &self,
    ) -> Result<AsyncSmtpTransport<Tokio1Executor>, lettre::transport::smtp::Error> {
        let credentials = Credentials::new(
            self.config.smtp.username.to_owned(),
            self.config.smtp.password.to_owned(),
        );

        let transport = AsyncSmtpTransport::<Tokio1Executor>::starttls_relay(
            &self.config.smtp.host.to_owned(),
        )?
            .port(self.config.smtp.port)
            .credentials(credentials)
            .build();

        Ok(transport)
    }

    fn render_template(&self, template_name: &str, template_data: &serde_json::Value) -> Result<String, handlebars::RenderError> {
        let mut handlebars = Handlebars::new();

        handlebars.register_template_file(template_name, &format!("/data/templates/email/{}.hbs", template_name))?;
        handlebars.register_template_file("styles", "/data/templates/email/styles.hbs")?;
        handlebars.register_template_file("base", "/data/templates/email/base.hbs")?;

        let content_template = handlebars.render(template_name, template_data)?;

        Ok(content_template)
    }

    async fn send_email(
        &self,
        subject: &str,
        template_name: &str,
        template_data: &serde_json::Value
    ) -> Result<(), Box<dyn std::error::Error>> {
        let html_template = self.render_template(template_name, template_data)?;
        let email = Message::builder()
            .to(
                format!("{} <{}>", self.user.name.as_str(), self.user.email.as_str())
                    .parse()
                    .unwrap(),
            )
            .reply_to(self.from.as_str().parse().unwrap())
            .from(self.from.as_str().parse().unwrap())
            .subject(subject)
            .header(ContentType::TEXT_HTML)
            .body(html_template)?;

        let transport = self.new_transport()?;
        
        transport.send(email).await?;

        Ok(())
    }

    pub async fn send_verification(&self, verification_url: &str, dev_version: bool) -> Result<(), Box<dyn std::error::Error>> {
        let template_data = TemplateData::new(
            "Account verification", 
            &self.user.name, 
            match dev_version { true => "Features may be unstable.", false => ""},
            verification_url,
            "(no signature provided)"
        );
        self.send_email(
            &template_data.subject,
            "verification",
             &template_data.to_json()
        ).await
    }
    
    pub async fn send_password_reset(&self, password_reset_url: &str) -> Result<(), Box<dyn std::error::Error>> {
        let template_data = TemplateData::new(
            "Password reset", 
            &self.user.name, 
            "",
            password_reset_url,
            "(no signature provided)"
        );
        self.send_email(
            &template_data.subject,
            "password",
             &template_data.to_json()
        ).await
    }
}
