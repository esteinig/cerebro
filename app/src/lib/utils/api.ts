import type { AuthRefreshResponseData, ErrorResponseData } from "./types";
import { env } from "$env/dynamic/public";
import { goto } from "$app/navigation";
import type {  ToastSettings } from "@skeletonlabs/skeleton";

export class ApiResponse {
    ok: boolean | undefined;
    status: number | undefined;
    statusText: string | undefined;
    json: any | null; 
    response: Response | null

    // If a response had the body already awaited
    // we need to return it separately - we mirror
    // some important field of the response in this
    // specific struct - this allows for refresh
    // and subsequent success/failure handling. We
    // also return the whole response in case it
    // is needed.
    constructor(response: Response | null, json: any | null, ok: boolean | null = null) {
        
        this.ok = ok === null ? response?.ok : ok;
        this.status = response?.status;
        this.statusText = response?.statusText;
        this.json = json;
        this.response = response;
    }
}
/**
 * Cerebro API 
 * 
 * API routes configured with public or private application configuration.
 * 
 * Public: uses base URL for public fetch requests from client-side browser to public API
 * Private: uses base URL for private fetch requests from server-side to Docker network API
 * 
 * @file lib/utils/config
 * @param {Routes} routes  Routes of configured Cerebro API
 */
export class CerebroApi {
    routes: Routes;

    constructor(url: string = env.PUBLIC_CEREBRO_API_URL) {
        this.routes = new Routes(url);
    }

   /**
     * Attempt to refresh the access token in client-side fetch calls
     * 
     * It may be that on this client-side fetch the access token expires
     * and a refresh token should be generated - however, the refresh logic
     * can only be implemented server-side (handle and handleFetch) - so we
     * need a helper function for each of these client-side defined requests
     * to protected routes, otherwise the page needs to be reloaded (which
     * is not ideal particularly for short lived access tokens) 
    */

    async fetchWithRefresh(url: string, init: RequestInit, refreshToken: string, toastStore: any | null = null, successMessage: string | null = "Success", timeout: number = 2000): Promise<ApiResponse> {
        

        try {
            let response: Response = await fetch(url, init);

            if (response.ok) { 
                let responseData: any = await response.json();
    
                if (toastStore !== null && successMessage !== null) {
                    toastStore.trigger({
                        message: successMessage,
                        timeout: timeout,
                        background: "variant-filled-primary",
                    } satisfies ToastSettings)
                };
    
                return new Promise<ApiResponse>((resolve) => { resolve(new ApiResponse(response, responseData)) })
            } else {
    
                // On error the response may contain a field 'refresh' which is returned by the
                // authentication middleware if the token could not be verified (due to expiration)
                let errorResponseData: ErrorResponseData;
    
                try {
                    errorResponseData = await response.json();
                } catch {
                    // Response may not have a body (Bad Request, Network Error)
                    if (toastStore !== null) {
                        toastStore.trigger({
                            message: `${response.status}: ${response.statusText}`,
                            timeout: timeout,
                            background: "variant-filled-tertiary",
                        } satisfies ToastSettings);
                    };
                    
                    return new Promise<ApiResponse>((resolve) => { resolve(new ApiResponse(response, {
                        status: `${response.status}`,
                        message: response.statusText
                    } as ErrorResponseData)) })
                }
    
                
                if (errorResponseData.refresh){
    
                    // In this case the access token cookie is expired in browser and not sent 
                    // with the request (but the refresh token is, if still valid)
    
                    // We set the refresh cookie directly - for some reason it is not sent
                    // with this request
                    let refreshResponse: Response = await fetch(this.routes.auth.refresh, {
                        method: 'GET',
                        mode: 'cors',
                        credentials: 'include'
                    } as RequestInit);
    
                    // Use the refreshed token directory
    
                    if (refreshResponse.ok){
    
                        // If the refresh is successful, rerun the original request and return response
                        let repeatResponse = await fetch(url, init);
                        let repeatResponseData: any = await repeatResponse.json();
    
                        if (toastStore !== null && successMessage !== null) {
                            toastStore.trigger({
                                message: successMessage,
                                timeout: timeout,
                                background: "variant-filled-primary",
                            } satisfies ToastSettings);
                        }
    
                        return new Promise<ApiResponse>((resolve) => { resolve(new ApiResponse(repeatResponse, repeatResponseData)) })
                    } else {
                        // If the refresh is not successful, user must complete login again
                        if (toastStore !== null) {
                            toastStore.trigger({
                                message: "Your session has expired. Please login again.",
                                background: "variant-filled-tertiary",
                                timeout: 3000
                            } satisfies ToastSettings);
                        }
                        setTimeout(() => goto("/login"), 3000);
                        return new Promise<ApiResponse>((resolve) => { resolve(new ApiResponse(response, errorResponseData)) })
                    }                
                } else {
                    // If a refresh indicator was not in the initial response,
                    // return the failed response
                    if (toastStore !== null) {
                        toastStore.trigger({
                            message: errorResponseData.message,
                            background: "variant-filled-tertiary",
                            timeout: timeout,
                        } satisfies ToastSettings);
                    }
    
                    return new Promise<ApiResponse>((resolve) => { resolve(new ApiResponse(response, errorResponseData)) })
                }
            }
        } catch {
            // Response may not have a body (Bad Request, Network Error)
            if (toastStore !== null) {
                toastStore.trigger({
                    message: `Failed to submit request`,
                    timeout: timeout,
                    background: "variant-filled-tertiary",
                } satisfies ToastSettings);
            };
            return new Promise<ApiResponse>((resolve) => { resolve(new ApiResponse(null, {
                status: '500',
                message: 'Critical failure to make a request'
            } as ErrorResponseData, false)) })
        }        
    }
}

interface AuthRoutes {
    login: string,
    logout: string,
    refresh: string,
    passwordEmail: string,
    verificationEmail: string,
    verificationCheck: string,
    passwordResetEmail: string,
    passwordResetCheck: string
}

interface UserRoutes {
    self: string,
    selfTeams: string,
    register: string,
    delete: string,
    update: string,
    get: string
}

interface TeamRoutes {
    register: string,
    delete: string,
    update: string,
    addUser: string,
    removeUser: string,
    get: string
}

interface CerebroRoutes {
    runs: string,
    workflows: string,
    samples: string,
    taxa: string,
    createReport: string,
    getReport: string,
    deleteReport: string,
    taxaEvidence: string,
    updateSampleDescription: string,
    deleteSamples: string,
    sampleOverview: string,
    getSampleSummary: string,
    searchSampleOverview: string,
    addSampleComment: string,
    deleteSampleComment: string,
    addPriorityTaxon: string,
    deletePriorityTaxon: string,
    updatePriorityTaxaDecision: string,
    updatePrioritTaxaDecisionComment: string
}

interface LogsRoutes {
    admin: string
}
interface FilesRoutes {
    getFiles: string,
    updateTags: string
}
interface PipelineRoutes {
    getPipelines: string
}

interface WatcherRoutes {
    getWatchers: string
}

interface StageRoutes {
    registerSamples: string
}


/**
 * Cerebro API Routes
 * 
 * @file lib/utils/config
 * @param {string} auth      auth routes
 * @param {string} users     user routes
 * @param {string} teams     team routes
 * @param {string} cerebro   cerebro routes
 */
export class Routes {
    auth: AuthRoutes;
    users: UserRoutes;
    teams: TeamRoutes;
    logs: LogsRoutes;
    files: FilesRoutes;
    pipelines: PipelineRoutes;
    watchers: WatcherRoutes;
    stage: StageRoutes;
    cerebro: CerebroRoutes;

    constructor(apiUrl: string) {

        let authRoute = `${apiUrl}/auth`;
        let userRoute = `${apiUrl}/users`;
        let teamRoute = `${apiUrl}/teams`;
        let logsRoute = `${apiUrl}/logs`;
        let filesRoute = `${apiUrl}/files`;
        let pipelinesRoute = `${apiUrl}/pipeline`;
        let watcherRoute = `${apiUrl}/watcher`;
        let stageRoute = `${apiUrl}/stage`;
        let cerebroRoute = `${apiUrl}/cerebro`;

        this.auth = {
            login: `${authRoute}/login`,
            logout: `${authRoute}/logout`,
            refresh: `${authRoute}/refresh`,
            verificationEmail: `${authRoute}/verification-email`,
            verificationCheck: `${authRoute}/verification-check`,
            passwordEmail: `${authRoute}/password-email`,
            passwordResetEmail: `${authRoute}/password-reset-email`,
            passwordResetCheck: `${authRoute}/password-reset-check`,
        };
        this.users = {
            self: `${userRoute}/self`,
            selfTeams: `${userRoute}/self/teams`,
            register: userRoute,
            delete: userRoute,
            update: userRoute,
            get: userRoute
        }
        this.teams = {
            register: teamRoute,
            delete: teamRoute,
            update: teamRoute,
            addUser: teamRoute,
            removeUser: teamRoute,
            get: teamRoute
        }
        this.cerebro = {
            runs: `${cerebroRoute}/runs`,
            workflows: `${cerebroRoute}/workflows`,
            samples: `${cerebroRoute}/samples`,
            taxa: `${cerebroRoute}/taxa`,
            createReport: `${cerebroRoute}/reports`,
            getReport: `${cerebroRoute}/reports`,
            deleteReport: `${cerebroRoute}/reports`,
            taxaEvidence: `${cerebroRoute}/taxa`,
            deleteSamples: `${cerebroRoute}/samples`,
            updateSampleDescription: `${cerebroRoute}/samples/description`,
            getSampleSummary: `${cerebroRoute}/samples/summary/qc`,
            sampleOverview: `${cerebroRoute}/samples/overview`,
            searchSampleOverview: `${cerebroRoute}/samples/overview`,
            addSampleComment: `${cerebroRoute}/samples/comment`,
            deleteSampleComment: `${cerebroRoute}/samples/comment`,
            addPriorityTaxon: `${cerebroRoute}/priority-taxa`,
            deletePriorityTaxon: `${cerebroRoute}/priority-taxa`,
            updatePriorityTaxaDecision: `${cerebroRoute}/priority-taxa/decision`,
            updatePrioritTaxaDecisionComment: `${cerebroRoute}/priority-taxa/decision/comment`,
        }
        this.logs = {
            admin: `${logsRoute}/admin`
        }
        this.files = {
            getFiles: `${filesRoute}`,
            updateTags: `${filesRoute}/tags`
        }
        this.pipelines = {
            getPipelines: `${pipelinesRoute}`
        }
        this.watchers = {
            getWatchers: `${watcherRoute}`
        }
        this.stage = {
            registerSamples: `${stageRoute}/register`
        }
    }
}


export default CerebroApi;