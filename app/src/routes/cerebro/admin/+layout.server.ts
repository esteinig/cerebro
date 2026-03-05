import type { User, ErrorResponseData, UsersResponseData, Team, TeamsResponseData, LogsResponseData, RequestLog } from "$lib/utils/types";
import { env as private_env } from "$env/dynamic/private";
import { error } from "@sveltejs/kit";
import CerebroApi from "$lib/utils/api";

let api = new CerebroApi(private_env.PRIVATE_CEREBRO_API_URL_DOCKER);

async function extractError(res: Response): Promise<string> {
    const contentType = res.headers.get("content-type") ?? "";
    const body = await res.text();

    if (contentType.includes("application/json")) {
        try {
            const parsed = JSON.parse(body);
            return parsed?.message ?? body;
        } catch {}
    }

    return body || `Request failed with status ${res.status}`;
}

const fetchUsers = async(fetch: Function, requestInit: RequestInit): Promise<Array<User>> => {
    
    let users: Array<User> = [];

    let usersResponse: Response = await fetch(
        `${api.routes.users.get}`, requestInit
    );
    if (usersResponse.ok) {
        let usersResponseData: UsersResponseData = await usersResponse.json();
        users = usersResponseData.data.users;
    } else {
        throw error(usersResponse.status, await extractError(usersResponse));
    }
    return users
}

const fetchTeams = async(fetch: Function, requestInit: RequestInit): Promise<Array<Team>> => {
    
    let teams: Array<Team> = [];

    let teamsResponse: Response = await fetch(
        `${api.routes.teams.get}`, requestInit
    );
    if (teamsResponse.ok) {
        let teamsResponseData: TeamsResponseData = await teamsResponse.json();
        teams = teamsResponseData.data.teams;
    } else {
        throw error(teamsResponse.status, await extractError(teamsResponse));
    }
    return teams
}


const fetchLogs = async(fetch: Function, requestInit: RequestInit): Promise<Array<RequestLog>> => {
    
    let logs: Array<RequestLog> = [];

    let logsResponse: Response = await fetch(
        `${api.routes.logs.admin}?limit=200`, requestInit
    );
    if (logsResponse.ok) {
        let logsResponseData: LogsResponseData = await logsResponse.json();
        logs = logsResponseData.data.logs;
    } else {
        throw error(logsResponse.status, await extractError(logsResponse));
    }
    return logs
}

const fetchCriticalLogs = async(fetch: Function, requestInit: RequestInit): Promise<Array<RequestLog>> => {
    
    let logs: Array<RequestLog> = [];

    let logsResponse: Response = await fetch(
        `${api.routes.logs.admin}?critical=true`, requestInit
    );
    
    if (logsResponse.ok) {
        let logsResponseData: LogsResponseData = await logsResponse.json();
        logs = logsResponseData.data.logs;
    } else {
        throw error(logsResponse.status, await extractError(logsResponse));
    }
    return logs
}

/** @type {import('./$types').PageLoad} */
export async function load({ fetch, depends }) {

    depends("admin:data")

    const fetchDataRequestInit: RequestInit = {
        method: 'GET',
        mode: 'cors',
        credentials: 'include'
    };

    return { users: fetchUsers(fetch, fetchDataRequestInit), teams: fetchTeams(fetch, fetchDataRequestInit), logs: fetchLogs(fetch, fetchDataRequestInit), criticalLogs: fetchCriticalLogs(fetch, fetchDataRequestInit) };
}