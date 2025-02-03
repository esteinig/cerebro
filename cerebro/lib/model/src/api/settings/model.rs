pub struct TeamSettings {
    pub report: ReportSettings
}

pub struct ReportSettings {
    pub templates: Vec<ReportTemplate>
}

pub struct ReportTemplate {
    pub template_id: String,
    pub template_name: String,
    pub contact_name: String,
    pub contact_email: String,
    pub legal_disclaimer: String,
    pub legal_disclosure: String,
    pub legal_liability: String,
    pub signatures: Vec<ReportSignature>
}

pub struct ReportSignature {
    pub name: String,
    pub position: String,
    pub institution: String
}

pub struct UserSettings {

}
