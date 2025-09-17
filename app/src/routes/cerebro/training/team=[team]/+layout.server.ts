import type { ErrorResponseData, Team, TrainingPrefetchOverview, TrainingResponse } from "$lib/utils/types";
import { error } from "@sveltejs/kit";
import type { PageServerLoad } from "../$types";
import { env as private_env } from "$env/dynamic/private";
import CerebroApi from "$lib/utils/api";

let api = new CerebroApi(private_env.PRIVATE_CEREBRO_API_URL_DOCKER);

const fetchOverview = async(fetch: Function, requestInit: RequestInit, team_param: string): Promise<TrainingPrefetchOverview[]> => {
    
    let training: TrainingPrefetchOverview[] = [];

    let trainingOverviewResponse: Response = await fetch(
        `${api.routes.training.getOverview}?team=${team_param}`, requestInit
    );
    if (trainingOverviewResponse.ok) {
        let trainingOverviewResponseData: TrainingResponse<TrainingPrefetchOverview[]> = await trainingOverviewResponse.json();
        training = trainingOverviewResponseData.data ?? [];
    } else {
        if (trainingOverviewResponse.status == 404) {  // No overview data found returns empty for page to render
            return training
        } else {
            let errorResponse: ErrorResponseData = await trainingOverviewResponse.json();
            throw error(trainingOverviewResponse.status, errorResponse.message)
        }
    }
    return training
}



export const load: PageServerLoad = async ({ params, locals, fetch, depends }) => { 

    depends("training:overview")

    const fetchDataRequestInit: RequestInit = {
        method: 'GET',
        mode: 'cors',
        credentials: 'include'
    };

    let team: Team;
    let currentUserTeams: Team[] = locals.teams;

    if (!currentUserTeams.length){
        throw error(404, `You are not a member of any teams - ${locals.admin ? "create a team for data upload first." : "please contact your administrator."}`)
    }

    // Special case on initial load of the page for default values
    if (params.team === "0") {
        team = currentUserTeams[0];
    } else {
        let selectedTeam = currentUserTeams.find(team => team.id === params.team || team.name === params.team);
        
        if (selectedTeam) { 
            team = selectedTeam 
        } else { 
            throw error(403, "You are not a member of this team") 
        }
    }
    
    let trainingOverview: TrainingPrefetchOverview[] = await fetchOverview(fetch, fetchDataRequestInit, team.name);

    return { 
        selectedTeam: team,
        trainingOverview: trainingOverview
    };

}