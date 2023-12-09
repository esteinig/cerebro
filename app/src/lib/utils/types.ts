
/**
 * =============
 * GENERAL TYPES
 * =============
*/

/**
 * Page source data returned from the load function in /routes/+layout.ts
 * 
 * This is available to all children routes and holds the SvelteKit fetch
 * function exposed by load.
 * 
 * @file lib/utils/types
 * @param {string} url - URL  of current route 
 * @param {string} commit - Latest git commit in repository (env var: PUBLIC_GIT_COMMIT) 
 */
export type allPageData = {
    url: string,
    commit: string,         
}

/**
 * Select option item type
 * 
 * @file lib/utils/types
 * @param {string} name - Label shown in select input
 * @param {number} value - Value associated with label
 */
export type SelectOption = {
    name: string, 
    value: string | number
}

/**
 * Unknown errors (exceptions) in the catch section of try-catch blocks
 * 
 * @file lib/utils/types
 * @param {name} status - Name of the exception 
 * @param {string} message - Status message of the exception
 */
export type UnknownException = {
    name: string,                       
    message: string
}

/**
 * Fetch request parameters usually used in client-side request
 * 
 * @file lib/utils/types
 * @param {string} url - Request URL
 * @param {RequestInit} requestInit - Request options
 */
export type FetchRequestParams = {
    url: string,                       
    requestInit: RequestInit
}


/**
 * Json response for expected endpoint failures
 * 
 * @file lib/utils/types
 * @param {string} status - Response status indicator
 * @param {string} message - Status response message from authentication server
 * @param {string} refresh - If token could not be verified an optional refresh field may be 
 */
export type ErrorResponseData = {
    status: string,                       
    message: string
    refresh?: boolean
}


/**
 * ============
 * NOTIFICATION
 * ============
*/

/**
 * Notification class to configure notification toasts
 * 
 * @file lib/utils/types
 * @param {string} visible - Notification toast is visible
 * @param {string} class - Notification toast class configuration
 * @param {string} message - Notification toast message text
 */
export type toastNotification = {
    visible: boolean,
    class: string,
    message: string
} 


/**
 * ===============
 * META-GP ROUTES
 * ===============
*/

/**
 * Page source data returned from the load function in /routes/meta-gp/+layout.server.ts
 *  
 * This data is available to all /meta-gp child routes - it contains the user data
 * which is returned from the server-sice layout data fetch, which serves as 
 * an authentication mechanism for the /meta-gp routes
 * 
 * @file lib/utils/types
 * @param {samples} url - Sam 
 */
export type mgpPageData = {
    userData: User,
    userTeams: Array<Team>
} & allPageData  // extends the layout source data containing fetch

/**
 * =========================
 * MGP-CEREBRO SAMPLE ROUTES
 * =========================
*/


/**
 * Page source data returned from the load function in /routes/meta-gp/samples/+page.server.ts
 *  
 * @file lib/utils/types
 * @param {Array<SampleOverviewData>} sampleOverview - Server-side requested sample overview data 
 * @param {Team} defaultDatabase - Server-side selected team  
 * @param {TeamDatabase} defaultDatabase - Server-side selected database 
 * @param {ProjectCollection} defaultProject - Server-side selected project 
 * @param {string} defaultNegativeTemplateControl - Server-side selected negative template control tag 
 * @param {number} defaultPageLimit - Server-side selected page limit
 */
export type CerebroSampleOverviewPageData = {
    // Server-side requested data
    sampleOverview: Array<SampleOverviewData>,
    // Server-side applied settings
    defaultTeam: Team,
    defaultDatabase: TeamDatabase,
    defaultProject: ProjectCollection,
    defaultNegativeTemplateControl: string,
    defaultPageLimit: number,
    workflowConfigs: Array<WorkflowConfig>,
    runConfigs: Array<RunConfig>,
    sampleGroups: Array<string>
} & mgpPageData  // extends the layout source data containing fetch


/**
 * Page source data returned from the load function in /routes/meta-gp/samples/+page.server.ts
 *  
 * @file lib/utils/types
 * @param {Array<Cerebro>} sampleCerebro - the sample specific data models (no taxa)
 * @param {Array<Cerebro>} controlCerebro - the sample associated run controls data models (no taxa)
 * @param {Array<WorkflowConfig>} sampleWorkflows - the sample associated workflow configurations
 * @param {WorkflowConfig} requestedWorkflow - the requested workflow configuration id
 */
export type CerebroSamplePageData = {
    // Server-side requested data
    sampleCerebro: Array<Cerebro>,
    controlCerebro: Array<Cerebro>,
    sampleWorkflows: Array<WorkflowConfig>,
    requestedWorkflow: string
} & mgpPageData  // extends the layout source data containing fetch


/**
 * Json response for successful sample overview response
 * 
 * @file lib/utils/types
 * @param {SampleOverviewResponseDataField} data - Sample overview response data
 */
export type SampleOverviewResponseData = {
    data: SampleOverviewResponseDataField
} & ErrorResponseData

/**
 * Json response for successful sample overview workflows response
 * 
 * @file lib/utils/types
 * @param {SampleOverviewResponseWorkflowsDataField} data - Sample overview workflows response data
 */
export type SampleOverviewWorkflowsResponseData = {
    data: SampleOverviewWorkflowsResponseDataField
} & ErrorResponseData


/**
 * Json response for successful sample overview response
 * 
 * @file lib/utils/types
 * @param {Array<SampleOverviewData>} sample_overview - Array of sample overview data items
 */
export type SampleOverviewResponseDataField = {
    sample_overview: Array<SampleOverviewData>
}

/**
 * Json response for successful sample workflows overview response
 * 
 * @file lib/utils/types
 * @param {Array<SampleOverviewWorkflowsData>} sample_overview - Array of sample overview workflows data items
 */
export type SampleOverviewWorkflowsResponseDataField = {
    sample_overview: Array<SampleOverviewWorkflowsData>
}

/**
 * Sample overview response model from the paginated aggregate pipeline at endpoint: /cerebro/samples/overview
 * 
 * @file lib/utils/types
 */
export type SampleOverviewData = {
    id: string,                           // sample identifier
    description: string[],           // sample descriptions (should be only one)
    latest_run: string,                   // latest workflow run date for this sample
    latest_workflow: string,              // latest workflow run date for this sample
    workflows: WorkflowConfig[],     // Workflow configurations for this sample
    samples: SampleConfig[],         // Sample configurations for this sample
    runs: RunConfig[],               // Run configurations for this sample
    tags: string[][],                // Library tags for all sample configs in this overview
    groups: string[],                // Groups for all sample configs in this overview
    types: string[],                 // Sampel types for all sample configs in this overview
    priority: PriorityTaxon[],
    reports: ReportEntry[]
}



/**
 * Quality control summary data from: /cerebro/samples/summary/qc
 * 
 * @file lib/utils/types
 */
export type QualityControlSummary = {
    id: string,                      // sample identifier
    model_id: string | null,         // model identifier
    sample_group: string | null,
    sample_type: string | null,
    run_id: string | null,
    run_date: string | null,
    workflow_name: string | null,
    workflow_id: string | null,
    workflow_date: string | null,
    total_reads: number,
    total_biomass: number | null, 
    deduplicated_reads: number,
    deduplicated_percent: number | null,
    ercc_constructs: number | null,
    ercc_reads: number | null,
    ercc_percent: number | null,
    ercc_input_mass: number | null,
    ercc_mass_per_read: number | null,
    adapter_trimmed_reads: number | null,
    adapter_trimmed_percent: number | null,
    low_complexity_reads: number | null,
    low_complexity_percent: number | null,
    mean_length_r1: number,
    mean_length_r2: number,
    qc_reads: number,
    qc_percent: number | null,
    qc_missing_bases_reads: number,
    qc_missing_bases_percent: number | null,
    qc_min_length_reads: number,
    qc_min_length_percent: number | null,
    qc_low_quality_reads: number,
    qc_low_quality_percent: number | null,
    q20_percent: number | null,
    q30_percent: number | null,
    host_reads: number | null,
    host_percent: number | null,
    host_biomass: number | null,
    other_reads: number | null,
    other_percent: number | null,
    other_biomass: number | null,
    dna_phage_id: string | null,
    dna_phage_reads: number | null,
    dna_phage_percent: string | null,
    dna_phage_coverage_percent: string | null,
    dna_phage_biomass: number | null,
    rna_phage_id: string | null,
    rna_phage_reads: number | null,
    rna_phage_percent: string | null,
    rna_phage_coverage_percent: string | null,
    rna_phage_biomass: number | null,
    seq_phage_id: string | null,
    seq_phage_reads: number | null,
    seq_phage_percent: string | null,
    seq_phage_coverage_percent: string | null,
    seq_phage_biomass: number | null,
    output_reads: number,
    output_percent: string | null,
    control_status_dna_extraction: boolean | null,
    control_status_rna_extraction:  boolean | null,
    control_status_library: boolean | null,
    control_status_sequencing: boolean | null,
}

/**
 * Quality control thresholds to indicate pass 
 * 
 * Will be set by user in team interface for preset storage in DB
 * 
 * @file lib/utils/types
 */
export type QualityControlThresholds = {
    input: Threshold[],
    deduplication: Threshold[],
    synthetic_controls: Threshold[],
    read_quality: Threshold[],
    host_depletion: Threshold[],
    phage_controls: Threshold[],
    output: Threshold[]
}

/**
 * Threshold value and logic with description to inform quality control status
 */
export type Threshold = {
    value: number,
    operator: string,
    description: string,
    symbol: string
}




/**
 * Sample overview response model from the match aggregate pipeline at endpoint: /cerebro/samples/overview/{id}
 * 
 * @file lib/utils/types
 * @param {}} id - Sample identifier
 */
export type SampleSummaryResponse = {
    data: SampleSummaryDataField
} & ErrorResponseData

/**
 * Sample quality control response model at endpoint: /cerebro/samples/summary/qc
 * 
 * @file lib/utils/types
 * @param {QualityControlSummary[]} summaries - Summaries
 * @param {string} csv - CSV string, empty string if not requested
 */
export type SampleSummaryDataField = {
    summary: QualityControlSummary[],
    csv: string
} & ErrorResponseData


/**
 * Sample overview response model from the match aggregate pipeline at endpoint: /cerebro/samples/overview/{id}
 * 
 * @file lib/utils/types
 * @param {string} id - Sample identifier
 * @param {Array<WorkflowConfig>} workflows - List of workflow configurations for this sample
 */
export type SampleOverviewWorkflowsData = {
    id: string,                    
    workflows: Array<WorkflowConfig>
}


/**
 * Sample delete schema to submit array of identifiers to delete: /cerebro/samples
 * 
 * @file lib/utils/types
 * @param {string[]} sample_id - Sample identifiers
 */
export type SampleDeleteSchema = {
    sample_id: string[],                    
}


/**
 * Sample summary schema to submit array of either sample identifiers or model identifiers to get quality control sumamries for: /cerebro/samples/summary/qc
 * 
 * @file lib/utils/types
 * @param {string[]} sample_id - Sample identifiers
 * @param {string[]} cerebro_id - Cerebro identifiers
 */
export type SampleSummarySchema = {
    sample_id: string[],         
    cerebro_id: string[]           
}


/**
 * Sample description schema to update a collections sample descriptions: /cerebro/samples/description
 * 
 * @file lib/utils/types
 * @param {string} description - Sample description
 */
export type SampleDescriptionSchema = {
    description: string,                
    sample_group: string,              
    sample_type: string          
}

/**
 * Json response for successful Cerebro model overview response (taxa removed)
 * 
 * @file lib/utils/types
 * @param {CerebroNotaxaDataField} data - Sample overview response data
 */
export type CerebroNotaxaResponse = {
    data: CerebroNotaxaDataField
} & ErrorResponseData



/**
 * Json response for successful Cerebro model response (taxa removed)
 * 
 * @file lib/utils/types
 * @param {Array<SampleOverviewData>} sample_overview - Array of sample overview data items
 */
export type CerebroNotaxaDataField = {
    cerebro: Array<Cerebro>  
}


/**
 * Json response for successful RunConfig model response 
 * 
 * @file lib/utils/types
 * @param {RunResponseDataField} data - Run config response data
 */
export type RunResponseData = {
    data: RunResponseDataField
} & ErrorResponseData




/**
 * Json response for successful RunConfig model response 
 * 
 * @file lib/utils/types
 * @param {Array<RunConfig>} runs - Array of unique run configs
 */
export type RunResponseDataField = {
    runs: Array<RunConfig>  
}



/**
 * Json response for successful WorkflowConfig model response 
 * 
 * @file lib/utils/types
 * @param {WorkflowResponseDataField} data - Workflow config response data
 */
export type WorkflowResponseData = {
    data: WorkflowResponseDataField
} & ErrorResponseData


/**
 * Json response for successful WorkflowConfig model response 
 * 
 * @file lib/utils/types
 * @param {Array<WorkflowConfig>} workflows - Array of unique workflow configs
 */
export type WorkflowResponseDataField = {
    workflows: Array<WorkflowConfig>  
}


/**
 * Json response for successful SampleConfig model response 
 * 
 * @file lib/utils/types
 * @param {SampleResponseData} data - Sample config response data
 */
export type SampleResponseData = {
    data: SampleResponseDataField
} & ErrorResponseData


/**
 * Json response for successful SampleConfig model response 
 * 
 * @file lib/utils/types
 * @param {Array<SampleConfig>} samples - Array of unique sample configs
 */
export type SampleResponseDataField = {
    samples: Array<SampleConfig>  
}

/*
=================
SAMPLE COMMENTS
=================
*/

export type SampleCommentSchema = {
    user_id: string,
    user_name: string,
    user_positions: string[],
    user_title: string | null,
    comment_text: string     
}

export type SampleComment = {
    comment_id: string,                          // assigned on upload
    user_id: string,
    user_name: string,
    user_title: string | null,
    user_positions: string[],
    comment_date: string,
    comment_text: string     
}


/**
 * ========================
 * AUTHENTICATION AND USERS
 * ========================
*/

/**
 * Page source data returned from the load function in /routes/meta-gp/settings/admin/+layout.server.ts
 * 
 * @file lib/utils/types
 * @param {Array<User>} users - the array of users returned from the /users endpoint (API) 
 * @param {Array<Team>} teams - the array of teams returned from the /teams endpoint (API) 
 */
export type adminSettingsPageData = {
    users: Array<User>,
    teams: Array<Team>
} & allPageData  // extends the layout source data containing fetch


/**
 * Schema to send one-time check token
 * 
 * @file lib/utils/types=
 * @param {string} access_token - Access token for one-time check routes
 */
export type AuthOneTimeCheckSchema = { 
    access_token: string
}

/**
 * Json response for successful token refresh server route
 * 
 * @file lib/utils/types
 * @param {string} status - Response status indicator
 * @param {string} access_token - Access token for local browser storage
 */
export type AuthRefreshResponseData = {
    status: string,    
    access_token: string
}

/**
 * Json response data for authentication endpoint failure
 * 
 * This extension of `ErrorResponse` contains an indicator whether to tryt refreshing
 * the access token in the server hook - this is usually when the token is expired
 * 
 * @file lib/utils/types
 * @param {boolean} refresh - Indicator whether to try refreshing the access token
 */
export type AuthRefreshErrorResponseData = {
    refresh: boolean,    
} & ErrorResponseData


/**
 * User login schema for login server endpoint
 * 
 * @file lib/utils/types
 * @param {string} email - User email address for login
 * @param {string} password - User password for login
 */
export type AuthLoginSchema = {
    email: string,                       
    password: string
}

/**
 * Json response for user login server route
 * 
 * Ages are returned for server-side cookie setting
 * 
 * @file lib/utils/types
 * @param {string} status - Response status indicator
 * @param {string} access_token - Access token
 * @param {string} refresh_token - Refresh token
 */
export type AuthLoginResponseData = {
    status: string,    
    access_token: string,
    refresh_token: string
}


/**
 * 
 * User authorization roles for controlled 
 * access to guarded endpoints
 * 
 * @file lib/utils/types
 */
export enum Role {
    User = "User",
    Admin  = "Admin",
    Bot = "Bot",
    Data = "Data",
    Report = "Report"
}

/**
 * Team configuration
 * 
 * @file lib/utils/types
 * @param {string} id - Team identifier (Uuid)
 * @param {string} name - Team name
 * @param {string} description - Short description of the team
 * @param {Array<TeamDatabase>} databases- Databases belonging to this team
 * @param {Array<string>} users - User identifiers belonging to this team
 */
export type Team = {
    id: string,
    name: string,
    description: string,
    databases: Array<TeamDatabase>,
    users: Array<string>
}

/**
 * MongoDB database setup for a `Team`
 * 
 * @file lib/utils/types
 * @param {string} id - Database identifier (Uuid)
 * @param {string} name - Database name
 * @param {string} database - Sanitized name of the database for connection
 * @param {string} description - Short description of the database
 * @param {Array<ProjectCollection>} projects - List of projects (collections) for this database
 */
export type TeamDatabase= {
    id: string,
    name: string,
    database: string,
    description: string,
    projects: Array<ProjectCollection>
}

/**
 * MongoDB collection setup for a `TeamDatabase`
 * 
 * @file lib/utils/types
 * @param {string} id - Project identifier (Uuid)
 * @param {string} name - Project name
 * @param {string} collection - Sanitized name of the collection for connection
 * @param {string} description - Short description of the database collection
 */
export type ProjectCollection = {
    id: string,
    name: string,
    collection: string,
    description: string,
    projects: Array<ProjectCollection>
}


/**
 * Team details for register/update inputs
 * 
 * @file lib/utils/types
 * @param {string} teamName - Team name
 * @param {string} teamLead - Unique user identifiers to register the team
 * @param {string} teamDescription - Short description of the team 
 * 
 */
export type TeamDetails = {
    teamName: string,
    teamLead: string,
    teamDescription: string,
}

/**
 * Team schema for registration server endpoint
 * 
 * @file lib/utils/types
 * @param {string} teamLead - Unique user identifiers to register the team
 * @param {string} teamName - Team name
 * @param {string} teamDescription - Short description of the team 
 * @param {string} databaseName - Team primary database name
 * @param {string} databaseMongoName - Sanitized name of the primary database in MongoDB
 * @param {string} databaseDescription - Short description of the primary database
 * @param {string} databaseName - Team primary database default project name
 * @param {string} databaseMongoName - Sanitized name of the primary database default project in MongoDB
 * @param {string} databaseDescription - Short description of the primary database defaualt project
 * 
 */
export type RegisterTeamSchema = {
    teamLead: string,
    teamName: string,
    teamDescription: string,
    databaseName: string,
    databaseMongoName: string,
    databaseDescription: string,
    projectName: string,
    projectMongoName: string,
    projectDescription: string,
}

/**
 * Team schema for update server endpoint
 * 
 * @file lib/utils/types
 * @param {string | null} teamName - Team name, null to not update
 * @param {string} teamDescription - Short description of the team
 * 
 */
export type UpdateTeamSchema = {
    teamName: string | null,
    teamDescription: string
}

/**
 * Status and data returned in successful team endpoint response (single team)
 * 
 * @file lib/utils/types
 * @param {string} status - Response status
 * @param {string} message - Response status message
 * @param {string} data - Response data with team
 */
export type TeamResponseData = {
    status: string,   
    message: string,                    
    data: TeamResponseDataField
}
/**
 * Team data returned in successful team endpoint response (single team)
 * 
 * @file lib/utils/types
 * @param {Team} team - Team data
 */
export type TeamResponseDataField = {                     
    team: Team
}


/**
 * User data filtered without password
 * 
 * @file lib/utils/types
 * @param {string} id - Identifier (Uuid)
 * @param {string} email - Email
 * @param {string} name - Full name 
 * @param {string | null} title- Title
 * @param {Array<string>} positions - Professional positions
 * @param {Array<Role>} roles - Resource roles
 * @param {string  | null} image - Profile image string
 * @param {boolean} verified - User is verified
 * @param {string} created - UTC ISO date 
 * @param {string} updated - UTC ISO date 
 */
export type User = {
    id: string,
    email: string,
    name: string,
    title: string | null,
    positions: Array<string>,
    roles: Array<Role>,
    image: string | null,
    verified: boolean,
    created: string,
    updated: string
}

/**
 * User register schema for registration server endpoint
 * 
 * @file lib/utils/types
 * @param {string} name - User first and last name
 * @param {string} email - User email address (must be unique)
 * @param {string} password - User password for registration
 * @param {boolean} verified - User is verified
 * @param {string | null} title - User preferred title
 * @param {Array<string>} positions - User professional positions
 * @param {Array<Role>} roles - Array of user access privileges
 */
export type RegisterUserSchema = {
    name: string,       
    email: string,       
    title: string | null,                
    password: string,
    verified: boolean,
    roles: Array<Role>,
    positions: Array<string>
}

/**
 * User data for update card assignment
 * 
 * @file lib/utils/types
 * @param {string} name - First and last name
 * @param {string} email - Email address (must be unique)
 * @param {string | null} password - Password for registration
 * @param {string  | null} title - Preferred title
 * @param {boolean} adminRole - Admin permissions
 * @param {boolean} dataRole - Data permissions
 * @param {boolean} reportRole - Report permissions
 * @param {boolean} botRole - Bot permissions
 * @param {string} position - Professional position
 * @param {boolean} position - Verified
 */
export type UpdateUserData = {
    name: string,       
    email: string,                   
    password: string | null,
    title: string | null,
    adminRole: boolean,
    dataRole: boolean,
    reportRole: boolean,
    botRole: boolean,
    position: string,
    verified: boolean
}

/**
 * User update schema for update server endpoint
 * 
 * @file lib/utils/types
 * @param {string} name - User first and last name
 * @param {string} email - User email address (must be unique)
 * @param {string | null} password - User password for registration
 * @param {string  | null} title - User preferred title
 * @param {Array<string>} positions - User professional positions
 * @param {boolean} verified - User is verified
 * @param {Array<Role>} roles - Array of access privileges
 */
export type UpdateUserScheme = {
    name: string,       
    email: string,                   
    password: string | null,
    title: string | null,
    roles: Array<Role>,
    verified: boolean,
    positions: Array<string>
}

/**
 * Status and data returned in successful user endpoint response (single user)
 * 
 * @file lib/utils/types
 * @param {string} data - Response data with user
 */
export type UserResponseData = {        
    data: UserResponseDataField
}

/**
 * User data returned in successful user endpoint response (single user)
 * 
 * @file lib/utils/types
 * @param {User} user - User data
 */
export type UserResponseDataField = {                     
    user: User
} & ErrorResponseData

/**
 * Status and data returned in successful user endpoint response (single user)
 * 
 * @file lib/utils/types
 * @param {string} data - Response data with user
 */
export type UserTeamsResponseData = {
    data: UserTeamsResponseDataField
} & ErrorResponseData

/**
 * User data returned in successful user endpoint response (single user)
 * 
 * @file lib/utils/types
 * @param {User} user - User data
 */
export type UserTeamsResponseDataField = {                     
    teams: Array<Team>
}

/**
 * Status and data returned in successful user endpoint response (multiple users)
 * 
 * @file lib/utils/types
 * @param {string} data - Response data containing users
 */
export type UsersResponseData = {           
    data: UsersResponseDataField
} & ErrorResponseData

/**
 * User data returned in successful multiple user endpoint response
 * 
 * @file lib/utils/types
 * @param {Array<User>} users - Array of users
 */
export type UsersResponseDataField = {                     
    users: Array<User>
}

/**
 * Status and data returned in successful user endpoint response (multiple teams)
 * 
 * @file lib/utils/types
 * @param {string} data - Response data containing teams
 */
export type TeamsResponseData = {              
    data: TeamsResponseDataField
} & ErrorResponseData

/**
 * Team data returned in successful multiple user endpoint response
 * 
 * @file lib/utils/types
 * @param {Array<Team>} teams - Array of teams
 */
export type TeamsResponseDataField = {                     
    teams: Array<Team>
}


/**
 * Status and data returned in successful logs endpoint response
 * 
 * @file lib/utils/types
 * @param {string} data - Response data containing teams
 */
export type LogsResponseData = {              
    data: LogsResponseDataField
} & ErrorResponseData

/**
 * Logs data returned in successful  logs endpoint response
 * 
 * @file lib/utils/types
 * @param {Array<Team>} teams - Array of teams
 */
export type LogsResponseDataField = {                     
    logs: Array<RequestLog>
}

/**
 * API request log entry
 * 
 * @file lib/utils/types
 * @param {string} id
 * @param {string} date 
 * @param {string} module
 * @param {string} action
 * @param {string} description 
 * @param {AccessDetails} access_details
 */
export type RequestLog = {
    id: string,
    date: string,
    module: string,
    action: string,
    description: string,
    access_details: AccessDetails,
}


/**
 * Details of a login access record
 * 
 * @file lib/utils/types
 * @param {RequestDetails} request_details
 * @param {string | null} user_id
 * @param {string | null} user_email 
 * @param {string | null} database_id
 * @param {string | null} project_id
 */
export type AccessDetails = {
    user_id: string | null,
    user_email: string | null,
    database_id: string | null,
    project_id: string | null,
    request_details: RequestDetails
}


/**
 * Request route details
 * 
 * @file lib/utils/types
 * @param {string} route
 * @param {string} method
 * @param {string[]} headers
 */
export type RequestDetails = {
    route: string,
    method: string,
    headers: string[]
}


/**
 * ===============
 * TAXONOMY MODELS
 * ===============
 */

export type Taxon = {
    taxid: string,
    rank: string,              // this is an enumeration in the `taxonomy` library, saved as string in model
    name: string,
    lineage: Array<string>,
    level: TaxonLevel,
    evidence: TaxonEvidence
}

export type TaxonLevel = {
    domain_taxid?: string,
    domain_name?: string,
    genus_taxid?: string,
    genus_name?: string,
    species_taxid?: string,
    species_name?: string
}

export type TaxonEvidence = {
    kmer: Kraken2UniqRecord[],
    alignment: VircovScanRemapRecord[],
    assembly: BlastLcaRecord[],
}

export type Kraken2UniqRecord = {
    id: string,
    db: string,
    tool: string,
    percent: string,
    reads: number,
    reads_direct: number,
    kmers: number,
    kmers_unique: number,
    tax_level: string,
    taxid: string,
    taxname: string,
    rpm: number
}


export type VircovRecord = {
    id: string,
    db: string,
    tool: string,
    taxid: string,
    reference: string,
    regions: number,
    reads: number,
    alignments: number,
    bases: number,
    length: number,
    coverage: number,
    description: string,
    tags: string,
    rpm: number
}

export type VircovScanRemapRecord = {
    id: string,
    db: string,
    tool: string,
    reference: string,
    reference_length: number,
    scan_regions: number,
    scan_reads: number,
    scan_alignments: number,
    scan_bases_covered: number,
    scan_coverage: number,
    remap_regions: number | null,
    remap_reads: number | null,
    remap_alignments: number | null,
    remap_bases_covered: number | null,
    remap_coverage: number | null,
    remap_mean_depth: number | null,
    consensus_length: number | null,
    consensus_missing: number | null,
    consensus_completeness: number | null,
    taxid: string | null,
    name: string | null,
    segment: string | null,
    reference_description: string,
    scan_rpm: number,
    remap_rpm: number,
}

export type BlastLcaRecord = {
    id: string,
    db: string,
    tool: string,
    length: number,    
    alignment: number,
    coverage: number,         
    identity: number,
    evalue: number,
    bitscore: number,
    taxid: string,
    title: string,
    reference: string,
    reference_length: number,
    read_coverage: number,           // if present in contig name e.g. from spades output
    bpm: number,
}

export type GenusOverview = {
    taxid: string,
    genus: string,                // used to later map back the tags
    domain: string | null,
    species: TaxonOverview[]
}

export type TaxonOverview = {
    taxid: string,
    name: string,                // used to later map back the tags
    domain: string | null,
    genus: string | null,
    rpm: number,                 // total rpm summed from k-mer and alignment evidence
    rpm_kmer: number,            // rpm from k-mer evidence
    rpm_alignment: number,       // rpm from alignment evidence
    contigs: number,             // total assembled and identified contig evidence
    contigs_bases: number,
    coverage: number,            // maximum coverage in the alignments evidence
    coverage_bases: number,      // maximum coverage bases of reference covered
    coverage_length: number,     // maximum coverage reference length
    kmer: boolean,
    alignment: boolean,
    assembly: boolean,
    names: Array<string>         // unique names of all evidence records (matching `Cerebro.name` field)
}


/**
 * Status and data returned in successful taxa overview endpoint response
 * 
 * @file lib/utils/types
 * @param {string} data - Response data containing overview taxa
 */
export type TaxaOverviewResponseData = {              
    data: TaxaOverviewResponseDataField
} & ErrorResponseData

/**
 * Taxon data returned in successful taxa overview response
 * 
 * @file lib/utils/types
 * @param {Array<TaxonOverview>} taxa - Array of taxa, summarized as overview models
 */

export type TaxaOverviewResponseDataField = {                     
    taxa: Array<TaxonOverview>
}

/**
 * =================
 * Frontend specific
 * =================
 */

export type PriorityTaxon = {
    id: string,
    user_name: string,
    user_id: string,
    date: string,
    comment: string,
    evidence_tags: Array<string>,  
    cerebro_identifiers: Array<string>,
    taxon_type: string,
    taxon_overview: TaxonOverview,
    filter_config: CerebroFilterConfig,
    decisions: Array<PriorityTaxonDecision>
};

export type PriorityTaxonSchema = {
    user_name: string,
    user_id: string,
    comment: string,
    evidence_tags: Array<string>,  
    cerebro_identifiers: Array<string>,
    taxon_type: PriorityTaxonType,
    taxon_overview: TaxonOverview,
    filter_config: CerebroFilterConfig,
    decisions: Array<PriorityTaxonDecision>
};

export type CandidateTaxonFormData = {
    userName: string,
    userIdentifier: string,
    comment: string,
    taxonType: PriorityTaxonType
}

/**
 * 
 * Decision type enumeration
 * 
 * @file lib/utils/types
 */
export enum DecisionType {
    Accept = "Accept",
    Reject  = "Reject"
}


/**
 * 
 * Candidate/priority taxon type enumeration
 * 
 * @file lib/utils/types
 */
export enum PriorityTaxonType {
    Unknown = "Unknown",
    Pathogen = "Pathogen",
    Contaminant = "Contaminant"
}

/**
 * 
 * Diagnostic status enumeration
 * 
 * @file lib/utils/types
 */
export enum DiagnosticStatus {
    Positive = "confirmed",
    Negative = "refuted",
    Unclear = "disputed",
    Unassigned = "unassigned"
}


type PriorityTaxonId = string;
type PriorityTaxonDecisionId = string;

export type PriorityTaxonDecision = { 
    id: PriorityTaxonDecisionId,
    date: string,
    user_id: string,
    user_name: string,
    decision: DecisionType,
    comment: string
}

export type PriorityTaxonDecisionSchema = { 
    id: PriorityTaxonId,
    decision: DecisionType,
    decisionComment: string,
    taxonName: string,      // For logging
    taxonType: string,      // For logging
    taxonTaxid: string,     // For logging
}

export type PriorityTaxonDecisionCommentSchema = { 
    priorityTaxonId: PriorityTaxonId,
    decisionId: PriorityTaxonDecisionId,
    decisionComment: string,
}


/**
 * Status and data returned in successful update of the priority taxon decision by a user
 * 
 * @file lib/utils/types
 * @param {string} data - Response data containing updated priority taxon decision
 */
export type PriorityTaxonDecisionResponseData = {              
    data: PriorityTaxonDecisionResponseDataField
} & ErrorResponseData

/**
 * Taxon data returned in successful taxa overview response
 * 
 * @file lib/utils/types
 * @param {PriorityTaxonDecision} decision - Updated priortiy taxon decision
 */

export type PriorityTaxonDecisionResponseDataField = {                     
    decision: PriorityTaxonDecision
}


export type CerebroFilterConfig = {  // ADD WORKFLOW ID
    domains: Array<string>,
    tags: Array<string>,
    kmer_min_reads: number,
    alignment_min_reads: number,
    alignment_min_bases: number,
    alignment_min_regions: number,
    alignment_min_coverage: number,
    alignment_min_ref_length: number,
    assembly_min_contig_length: number,
    assembly_min_contig_identity: number
    assembly_min_contig_coverage: number
}

// Client-side selection of taxa

export type ClientFilterConfig = {
    domains: Array<string | null>,
    genera: Array<string | null>,
    species: Array<string | null>,
    modules: ClientFilterModules,
    minimum: ClientFilterMinimum
}

export type ClientFilterModules = {
    alignment: boolean,
    kmer: boolean,
    assembly: boolean
}
export type ClientFilterMinimum = {
    rpm: number,
    rpm_kmer: number,
    rpm_alignment: number,
    contigs: number,
    bases: number
}



/*
===========================
DATABASE CEREBRO BASE MODEL
===========================
*/

export type Cerebro = {
    schema_version: string,             // ther schema version of this model

    id: string,                         // the unique identifier of this model in the database
    name: string,                       // the identifier for the workflow sample used in the pipeline

    run: RunConfig,                     // the configuration of the sequence run this sample was sequenced on
    sample: SampleConfig,               // the configuration of the biological sample this sample was processed with
    workflow: WorkflowConfig,           // the configuration of the workflow run this sample was processed with 

    taxa?: Map<string, Taxon>,           // the dictionary of taxonomic identifiers and taxon data from the parsed workflow sample (legacy dictionary, could be simple list)
    quality: QualityControlModule,       // the quality control data from the parsed workflow sample

}

export type RunConfig = {               // parsed from the `sample_sheet.csv` of the nextflow
    id: string,                         // sequence run identifier set in sample sheet
    date: string                        // sequence run date set in the sample sheet
}

export type SampleConfig = {
    id: string,                         // biological sample identifier parse from input file name
    tags: Array<string>,                // library tags parsed from input file name
    description: string,                // workflow description
    sample_group: string,               // sample group defined in `sample_sheet.csv`
    sample_type: string,                // sample type defined in `sample_sheet.csv`
    priority: Array<PriorityTaxon>,      // priority taxa set from the user interface
    reports: Array<ReportEntry>,
    comments: Array<SampleComment>
}

export type ReportEntry = {
    id: string,
    date: string,
    user_id: string,
    user_name: string,
    negative: boolean,
    organism: string,
    review_date: string,
    report_text: string | null,
    report_pdf: boolean | null
}

export type WorkflowConfig = {          // parsed from `config.json` nextflow output
    id: string,                         // unique identifier of the workflow run 
    name: string,                       // mnenomic name of the workflow
    pipeline: string,
    version: string,                    // version of the workflow 
    started: string,                    // datetime of workflow run start 
    completed: string,                  // datetime of workflow run end 
    description: string,                // description of workflow
    params: WorkflowParams              // workflow parameters from 
}

export type WorkflowParams = {
    production: boolean                 // workflow run in production mode
    qc: WorkflowParamsQc
}

export type WorkflowParamsQc = {
    deduplication: WorkflowParamsQcDeduplication,
    reads: WorkflowParamsQcReads,
    controls: WorkflowParamsQcControls,
    host: WorkflowParamsQcHost
}

export type WorkflowParamsQcDeduplication = {
    enabled: boolean,
    method: string,
    umi_tools: WorkflowParamsQcUmiTools
}

export type WorkflowParamsQcUmiTools = {
    seed: number,
    reference: string
}

export type WorkflowParamsQcReads = {
    fastp: WorkflowParamsQcFastp,
}

export type WorkflowParamsQcFastp  = {
    enabled: boolean,
    min_read_length: number | null,
    cut_tail_quality: number | null,
    complexity_threshold: number | null,
    adapter_auto_detect: boolean,
    adapter_file: string | null,
    adapter_seq_1: string | null,
    adapter_seq_2: string | null,
    trim_poly_g: number | null,
}

export type WorkflowParamsQcControls  = {
    ercc: WorkflowParamsQcErcc
    phage: WorkflowParamsQcPhage
}

export type WorkflowParamsQcErcc  = {
    enabled: boolean,
    fasta: string | null
}
export type WorkflowParamsQcPhage  = {
    enabled: boolean,
    fasta: string | null,
    identifiers: WorkflowParamsQcPhageIdentifiers
}

export type WorkflowParamsQcPhageIdentifiers  = {
    dna_extraction: string,
    rna_extraction: string,
    sequencing: string
}

export type WorkflowParamsQcHost = {
    depletion: WorkflowParamsQcHostDepletion
}

export type WorkflowParamsQcHostDepletion = {
    enabled: boolean,
    databases: string,
    references: string,
    taxa: string,
    direct: string,
    min_cov: number,
    min_len: number,
    min_mapq: number
}

/*
=================
DATABASE QC MODEL
=================
*/

export type QualityControlModule = {
    id: string,
    fastp: Fastp,
    ercc: Ercc | null,
    host_background: Scrubby | null,
    virus_background: Scrubby | null
}

export type Fastp = {
    summary: FastpSummary,
    filtering_result: FastpFilter,
    duplication: FastpDuplication,
    adapter_cutting: FastpAdapters
}

export type FastpSummary = {
    fastp_version: string,
    before_filtering: FastpReadSummary,
    after_filtering: FastpReadSummary
}

export type FastpReadSummary  = {
    total_reads: number,
    total_bases: number,
    q20_rate: number,
    q30_rate: number,
    read1_mean_length: number,
    read2_mean_length: number,
    gc_content: number
}

export type FastpFilter = {
    passed_filter_reads: number,
    low_quality_reads: number,
    low_complexity_reads: number,
    too_many_N_reads: number,
    too_short_reads: number,
    too_long_reads: number
}

export type FastpDuplication = {
    rate: number
}

export type FastpAdapters = {
    adapter_trimmed_reads: number,
    adapter_trimmed_bases: number,
    read1_adapter_sequence: string,
    read2_adapter_sequence: string
}


export type Ercc = {
    constructs: number,
    constructs_aligned: number,
    reads: number,
    reads_aligned: number,
    depletion: Scrubby,
    records: Array<VircovRecord>
}

export type Scrubby = {
    version: string,
    schema_version: string,
    summary: ScrubbySummary,
    settings: ScrubbySettings,
    pipeline: Array<ScrubbyReferenceSummary>
}

export type ScrubbySummary = {
    total: number,
    depleted: number,
    extracted: number
}

export type ScrubbySettings = {
    kraken_taxa: Array<string>,
    kraken_taxa_direct: Array<string>,
    min_len: number,
    min_cov: number,
    min_mapq: number,
    extract: boolean
}

export type ScrubbyReferenceSummary = {
    index: number,
    tool?: string,
    name: string,
    path: string,
    total: number,
    depleted: number,
    extracted: number,
    files: Array<ScrubbyFileSummary>
}

export type ScrubbyFileSummary = {
    total: number,
    depleted: number,
    extracted: number,
    input_file: string,
    output_file: string,
}



/**
 * Comment bubble
 * 
 * @file lib/utils/types
 * @param {number} id
 * @param {boolean} host
 * @param {number} avatar
 * @param {string} name
 * @param {string} timestamp
 * @param {string} message
 * @param {string} color
 */
export type CommentBubble = {
    id: string,
    host: boolean,
    name: string,
    timestamp: string,
    message: string,
    color: string,
    title: string | null,
    position: string
}


/**
 * Report schema
 * 
 * @file lib/utils/types
 */
export type ReportSchema = {
    ids: string[],
    user_id: string,
    user_name: string,
    sample_id: string,
    run_id: string,
    negative: boolean,
    workflow: WorkflowConfig,
    patient_header: PatientHeaderSchema,
    patient_result: PatientResultSchema,
    laboratory_comments: string,
    bioinformatics_comments: string,
    priority_taxon: PriorityTaxon | null,
    priority_taxon_negative_control: boolean | null,
    priority_taxon_other_evidence: string | null,
}

export type PatientHeaderSchema = {
    patient_name: string,
    patient_urn: string,
    patient_dob: string,
    requested_doctor: string,
    hospital_site: string,
    laboratory_number: string,
    specimen_id: string,
    date_collected: string,
    date_received: string,
    specimen_type: string,
    reporting_laboratory: string,
    reporting_date: string,
    reporting_location: string,
}

export type PatientResultSchema = {
    review_date: string,
    negative: boolean,
    organism: string,
    contact: string,
    comments: string,
    actions: string,
}

export type TimelineField = {
    
}