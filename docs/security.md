This documentation describes the configuration and considerations around data privacy and security in `Cerebro`.

We take particular care in explaining how different models and configurations impact on data security of the deployment. This may seem overly cautious - after all we are quite used to things just working silently in the background.

However, we think this topic is extremely important and has not been sufficiently considered for the translation of research-focused metagenomics diagnostic methods into clinical production applications. This is particularly due to issues around [depletion and identifiability]() of host background sequences.
    
Regulatory and accreditation authorities may not always have the experience to consider the nuances of modern genomics implementations, infrastructure and deployments. It is the responsibility of the administrators of such systems to ensure that data is kept secure and conforms to local ethics and legal frameworks. 
    
We hope that this guide provides sufficient background information and alleviates uncertainties about how potentially sensitive data are handled in `Cerebro`. 

!!! info
    
    We favour a metagenomic diagnostic model that does not involve upload of sensitive data to servers controlled by third-parties with opaque data-privacy implementations or questionable motivations regarding data collection and storage.  We also want to give users the freedom to decide how and with what level of security they can deploy the `Cerebro` stack on their own infrastructure.

### TLDR

Sequence and patient data are de-coupled from the database and are not linked to the analytical results. However, there are user-determined points of caution and less restrictive configurations that can be activated by the stack administrator. These allow for data linkage within a secure network, for example to generate clinical reports with patient details.

We provide multiple configurations and let administrators decide what level of security they wish to implement for their deployment. It is the responsibility of users to evaluate and comply with local regulations ([Disclaimer]() and [Liability]()).


## Data Security

### Concepts

`Cerebro` de-couples the pipeline execution stage from storing, accessing and reporting analysis outputs. Sequence data can therefore be processed and analysed on a server within the network where the data was produced, for example within the secure network of the public health laboratory that conducts the sequencing operations.

Results of quality control, taxonomic classifications and accessory analyses are aggregated into serializable database models. Besides storing models in the database, they can be output as `JSON` files and are generally safe (see below) to share outside the secure network. Results contained in the models can be examined and further analysed through the user interface or command-line client.

### Implications

**No sequence data or patient details are stored in the database**. Sequence data or patient details are not linked to results and are not accessible with default configurations of `Cerebro`.


!!! warning

    There are three potential vulnerabilities for accidental data linkage that should be considered:

    1. Library (sample) identifiers should be de-identified and should not replicate identifiers used for patient details stored elsewhere. Sample identifiers and file paths that include the sample identifiers (such as input read files) are stored in the database.
    2. Additional annotations of samples can be added through the sample sheet or user interface and are stored in the database models - careful users should consider what meta-data about the sample is acceptable to be included.
    3. Comments by team members can be submitted through the user-interface and are stored in the database models. If team members discuss results in the context of patient details these could be linked to diagnostic results in the application. 
    
    Depending on where the stack is deployed and who administers the stack, annotations and comment functions may constitute an unacceptable risk or may not comply with local regulatory environments. Annotations and comments can therefore be [deactivated]() application-wide by the stack administrator during deployment.

#### Sequence Data

Processed sequence data is generally only accessible on the server where the pipeline was executed. 

#### Report Headers

Patient details are usually stored in secure and regulated infrastructure outside the scope of this application (e.g. REDCap in Australia). However, for the final clinical report, there is a necessity to include sensitive information in the header section of the clinical report returned to the clinician. 


We provide three options for handling the inclusion of such information depending on the desired level of security and compliance with local regulatory environments:

1. Reports can be generated without header information. This forces users to fill-in patient details in the resulting PDF.
2. Header information should only be provided through the application, if both the server and the interface are deployed on a local machine within a secure network, so that header data is never transmitted on the network. This assumes that the local machine and network are secure.
3. If the application server is deployed on a different machine within the secure network or on a remote server, any header information submitted through the application should be compiled with the web-assembly option. This forces compilation of the reports in the browser, so that patient information never leaves the local machine. This assumes that the local machine and browser are secure.

!!! warning

    We provide these options for users to select their desired level of security, as the purpose of deployment and the environment in which the application is run may be highly variable between sites and regulatory environments. 
    
    Depending on where the stack is deployed and who administers the stack, submitting any patient information to the application may not comply with legal frameworks and constitute an unacceptable risk. Reporting headers can therefore be [deactivated]() globally by the stack administrator during deployment, which forces users to fill-in the headers in the generated PDF.


    
    



## Disclaimer