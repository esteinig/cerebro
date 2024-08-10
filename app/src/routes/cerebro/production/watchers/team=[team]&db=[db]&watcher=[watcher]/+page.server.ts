import { type ErrorResponseData, type Team, type TeamDatabase, type SeaweedFile, type FileResponseData, type ProjectCollection, type PipelineResponseData, type ProductionPipeline, type ProductionWatcher, type WatcherResponseData, parseWatcherFormat, parsePipeline } from "$lib/utils/types";
import { env as private_env } from "$env/dynamic/private";
import { error } from "@sveltejs/kit";
import CerebroApi from "$lib/utils/api";
import type { PageServerLoad } from "./$types";

let api = new CerebroApi(private_env.PRIVATE_CEREBRO_API_URL_DOCKER);

const fetchFiles = async(fetch: Function, requestInit: RequestInit, db_param: string, watcher_param: string | undefined): Promise<SeaweedFile[]> => {
    
    let files: SeaweedFile[] = [];

    if (watcher_param === undefined) {
        return files
    }

    let filesResponse: Response = await fetch(
        `${api.routes.files.getFiles}?db=${db_param}&watcher_id=${watcher_param}&page=0&limit=1000`, requestInit
    );
    if (filesResponse.ok) {
        let filesResponseData: FileResponseData = await filesResponse.json();
        files = filesResponseData.data;
    } else {
        if (filesResponse.status == 404) {  // No files found returns empty for page to render
            return files
        } else {
            let errorResponse: ErrorResponseData = await filesResponse.json();
            throw error(filesResponse.status, errorResponse.message)
        }
    }
    return files
}


const fetchPipelines = async(fetch: Function, requestInit: RequestInit, db_param: string): Promise<ProductionPipeline[]> => {
    
    let pipelines: ProductionPipeline[] = [];

    let pipelinesResponse: Response = await fetch(
        `${api.routes.pipelines.getPipelines}?db=${db_param}`, requestInit
    );
    if (pipelinesResponse.ok) {
        let pipelinesResponseData: PipelineResponseData = await pipelinesResponse.json();

        pipelines = pipelinesResponseData.data.map(pipeline => ({
            ...pipeline,
            format: parsePipeline(pipeline.pipeline)
        }));
    } else {
        if (pipelinesResponse.status == 404) {  // No pipelines found returns empty for page to render
            return pipelines
        } else {
            let errorResponse: ErrorResponseData = await pipelinesResponse.json();
            throw error(pipelinesResponse.status, errorResponse.message)
        }
    }
    return pipelines
}

const fetchWatchers = async(fetch: Function, requestInit: RequestInit, db_param: string): Promise<ProductionWatcher[]> => {
    
    let watchers: ProductionWatcher[] = [];

    let watchersResponse: Response = await fetch(
        `${api.routes.watchers.getWatchers}?db=${db_param}`, requestInit
    );
    if (watchersResponse.ok) {
        let watchersResponseData: WatcherResponseData = await watchersResponse.json();

        watchers = watchersResponseData.data.map(watcher => ({
            ...watcher,
            format: parseWatcherFormat(watcher.format)
        }));
    } else {
        if (watchersResponse.status == 404) {  // No watchers found returns empty for page to render
            return watchers
        } else {
            let errorResponse: ErrorResponseData = await watchersResponse.json();
            throw error(watchersResponse.status, errorResponse.message)
        }
    }
    return watchers
}


export const load: PageServerLoad = async ({ params, locals, fetch, depends }) => {

    depends("watchers:data")

    const fetchDataRequestInit: RequestInit = {
        method: 'GET',
        mode: 'cors',
        credentials: 'include'
    };

    let currentUserTeams: Team[] = locals.teams;

    // Pick the first team and database from the user 
    // return an error if they are not part of a team
    // (which always has one default database)

    if (!currentUserTeams.length){
        throw error(404, `You are not a member of any teams - ${locals.admin ? "create a team for data upload first." : "please contact your administrator."}`)
    }

    let defaultTeam: Team = currentUserTeams[0];
    
    if (!defaultTeam.databases.length){
        throw error(404, `No database has been created for team (${defaultTeam.name})`)
    }

    let defaultDatabase: TeamDatabase = defaultTeam.databases[0];

    if (!defaultDatabase.projects.length){
        throw error(404, `No database projects have been created for team (${defaultTeam.name}) and database (${defaultDatabase.name})`)
    }
    let defaultProject: ProjectCollection = defaultDatabase.projects[0];
        
    const [registeredPipelines, registeredWatchers] = await Promise.all([
        fetchPipelines(fetch, fetchDataRequestInit, params.db),
        fetchWatchers(fetch, fetchDataRequestInit, params.db),
    ]);

    let defaultWatcher: ProductionWatcher | undefined = params.watcher === "0" ? registeredWatchers[0] : registeredWatchers.find(watcher => watcher.id === params.watcher);

    // Page initialisation without a selected watcher identifier since
    // we are getting them in this function - use the default watcher

    let files: SeaweedFile[] = await fetchFiles(
        fetch, 
        fetchDataRequestInit, 
        params.db, 
        params.watcher === "0" ? defaultWatcher?.id : params.watcher
    );

    return { 
        files: files,
        registeredPipelines: registeredPipelines,
        registeredWatchers: registeredWatchers,
        defaultTeam: defaultTeam,
        defaultDatabase: defaultDatabase,
        defaultProject: defaultProject,
        defaultWatcher: defaultWatcher
    };

}