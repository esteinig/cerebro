import type { ErrorResponseData, Team, TeamDatabase, SeaweedFile, FileResponseData } from "$lib/utils/types";
import { env as private_env } from "$env/dynamic/private";
import { error } from "@sveltejs/kit";
import CerebroApi from "$lib/utils/api";
import type { PageServerLoad } from "./$types";



/**
 * Group SeaweedFile objects by a specified key (watcher_location, run_id, or date).
 * Adjusted to handle nullable run_id by excluding those files from run_id grouping
 * or placing them under a specific key.
 *
 * @param files Array of SeaweedFile objects to be grouped.
 * @param groupByKey The key by which the files should be grouped ('watcher_location', 'run_id', or 'date').
 * @param dateRange (Optional) Specifies the date range for grouping if 'date' is chosen as groupByKey. 
 *                  Format: { start: 'YYYY-MM-DD', end: 'YYYY-MM-DD' }
 * @returns An object with keys as the group by value and values as arrays of SeaweedFiles belonging to that group.
 *
 * @example
 * // Example usage for grouping by watcher location
 * const files: SeaweedFile[] = [...];
 * const groupedByLocation = groupSeaweedFiles(files, 'watcher_location');
 *
 * // Example usage for grouping by run_id, excluding null run_ids
 * const groupedByRunId = groupSeaweedFiles(files, 'run_id');
 *
 * // Example usage for grouping by date with a specified range
 * const groupedByDateRange = groupSeaweedFiles(files, 'date', { start: '2022-01-01', end: '2022-01-31' });
 */
function groupSeaweedFiles(files: SeaweedFile[], groupByKey: 'watcher_location' | 'run_id' | 'date', dateRange?: { start: string, end: string }): Record<string, SeaweedFile[]> {
    const grouped: Record<string, SeaweedFile[]> = {};

    files.forEach(file => {
        let key: string | null = null;

        switch (groupByKey) {
            case 'watcher_location':
                key = file.watcher.location;
                break;
            case 'run_id':
                // Handle nullable run_id by setting a default key or excluding
                key = file.run_id ? file.run_id : 'no_run_id'; // Or `null` to exclude
                break;
            case 'date':
                if (dateRange && file.date >= dateRange.start && file.date <= dateRange.end) {
                    key = file.date;
                } else if (!dateRange) {
                    key = file.date;
                } // Files outside the specified range are excluded
                break;
            default:
                throw new Error(`Unsupported groupByKey: ${groupByKey}`);
        }

        if (key) {
            grouped[key] = grouped[key] || [];
            grouped[key].push(file);
        }
    });

    return grouped;
}

let api = new CerebroApi(private_env.PRIVATE_CEREBRO_API_URL_DOCKER);

const fetchFiles = async(fetch: Function, requestInit: RequestInit, db_param: string): Promise<SeaweedFile[]> => {
    
    let files: SeaweedFile[] = [];

    let filesResponse: Response = await fetch(
        `${api.routes.files.getFiles}?db=${db_param}&page=0&limit=1000`, requestInit
    );
    if (filesResponse.ok) {
        let filesResponseData: FileResponseData = await filesResponse.json();
        files = filesResponseData.data;
    } else {
        let errorResponse: ErrorResponseData = await filesResponse.json();
        throw error(filesResponse.status, errorResponse.message)
    }
    return files
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

    let team: Team = currentUserTeams[0];

    if (!team.databases.length){
        throw error(404, `No database has been created for this team (${team.name})`)
    }

    let teamDatabase: TeamDatabase = team.databases[0];

    let files = await fetchFiles(fetch, fetchDataRequestInit, params.db);

    let filesWatcherLocation = groupSeaweedFiles(files, "watcher_location");
    let defaultWatcherLocation: string = "";
    let keys = Object.keys(filesWatcherLocation);
    if (keys.length) {
        defaultWatcherLocation = keys[0];
    }

    return { 
        files: filesWatcherLocation,
        defaultTeam: team,
        defaultDatabase: teamDatabase,
        defaultWatcherLocation: defaultWatcherLocation,
     };
}