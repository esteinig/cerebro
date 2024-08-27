import { type ErrorResponseData, type Team, type TeamDatabase, type SeaweedFile, type FileResponseData, type ProjectCollection, type TowerResponseData, type ProductionTower, type ProductionWatcher, type WatcherResponseData, parseWatcherFormat, parsePipeline } from "$lib/utils/types";
import { env as private_env } from "$env/dynamic/private";
import { error } from "@sveltejs/kit";
import CerebroApi from "$lib/utils/api";
import type { PageServerLoad } from "./$types";

let api = new CerebroApi(private_env.PRIVATE_CEREBRO_API_URL_DOCKER);

const fetchFiles = async(fetch: Function, requestInit: RequestInit, team_param: string, watcher_param: string | undefined): Promise<SeaweedFile[]> => {
    
    let files: SeaweedFile[] = [];

    if (watcher_param === undefined) {
        return files
    }

    let filesResponse: Response = await fetch(
        `${api.routes.files.getFiles}?team=${team_param}&watcher_id=${watcher_param}&page=0&limit=1000`, requestInit
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


const fetchTowers = async(fetch: Function, requestInit: RequestInit, team_param: string): Promise<ProductionTower[]> => {
    
    let towers: ProductionTower[] = [];

    let towersResponse: Response = await fetch(
        `${api.routes.towers.getPipelines}?team=${team_param}`, requestInit
    );
    if (towersResponse.ok) {
        let towersResponseData: TowerResponseData = await towersResponse.json();

        towers = towersResponseData.data.map(tower => ({
            ...tower,
            format: tower.pipelines.map(pipeline => parsePipeline(pipeline))
        }));
    } else {
        if (towersResponse.status == 404) {  // No towers found returns empty for page to render
            return towers
        } else {
            let errorResponse: ErrorResponseData = await towersResponse.json();
            throw error(towersResponse.status, errorResponse.message)
        }
    }
    return towers
}

const fetchWatchers = async(fetch: Function, requestInit: RequestInit, team_param: string): Promise<ProductionWatcher[]> => {
    
    let watchers: ProductionWatcher[] = [];

    let watchersResponse: Response = await fetch(
        `${api.routes.watchers.getWatchers}?team=${team_param}`, requestInit
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

    let selectedTeam: Team = currentUserTeams.find(team => team.id === params.team) || currentUserTeams[0];
    
    if (!selectedTeam.databases.length){
        throw error(404, `No database has been created for team (${selectedTeam.name})`)
    }

    let defaultDatabase: TeamDatabase = selectedTeam.databases[0];

    if (!defaultDatabase.projects.length){
        throw error(404, `No database projects have been created for team (${selectedTeam.name}) and database (${defaultDatabase.name})`)
    }
    let defaultProject: ProjectCollection = defaultDatabase.projects[0];
        
    const [registeredTowers, registeredWatchers] = await Promise.all([
        fetchTowers(fetch, fetchDataRequestInit, params.team),
        fetchWatchers(fetch, fetchDataRequestInit, params.team),
    ]);

    let defaultWatcher: ProductionWatcher | undefined = params.watcher === "0" ? registeredWatchers[0] : registeredWatchers.find(watcher => watcher.id === params.watcher);

    // Page initialisation without a selected watcher identifier since
    // we are getting them in this function - use the default watcher

    let files: SeaweedFile[] = await fetchFiles(
        fetch, 
        fetchDataRequestInit, 
        params.team,
        params.watcher === "0" ? defaultWatcher?.id : params.watcher
    );

    return { 
        files,
        registeredTowers,
        registeredWatchers,
        selectedTeam,
        defaultDatabase,
        defaultProject,
        defaultWatcher
    };

}