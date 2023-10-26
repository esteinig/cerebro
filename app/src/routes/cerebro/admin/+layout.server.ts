import type { User, ErrorResponseData, UsersResponseData, Team, TeamsResponseData, LogsResponseData, RequestLog } from "$lib/utils/types";
import { env as private_env } from "$env/dynamic/private";
import { error } from "@sveltejs/kit";
import CerebroApi from "$lib/utils/api";

let api = new CerebroApi(private_env.PRIVATE_CEREBRO_API_URL_DOCKER);

const fetchUsers = async(fetch: Function, requestInit: RequestInit): Promise<Array<User>> => {
    
    let users: Array<User> = [];

    let usersResponse: Response = await fetch(
        `${api.routes.users.get}`, requestInit
    );
    if (usersResponse.ok) {
        let usersResponseData: UsersResponseData = await usersResponse.json();
        users = usersResponseData.data.users;
    } else {
        let errorResponse: ErrorResponseData = await usersResponse.json();
        throw error(usersResponse.status, errorResponse.message)
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
        let errorResponse: ErrorResponseData = await teamsResponse.json();
        throw error(teamsResponse.status, errorResponse.message)
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
        let errorResponse: ErrorResponseData = await logsResponse.json();
        throw error(logsResponse.status, errorResponse.message)
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
        let errorResponse: ErrorResponseData = await logsResponse.json();
        throw error(logsResponse.status, errorResponse.message)
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