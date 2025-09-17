import { env as private_env } from "$env/dynamic/private";
import CerebroApi from "$lib/utils/api";
import type { CreateTrainingSession, ErrorResponseData, Team, TrainingPrefetchData, TrainingRecord, TrainingResponse } from "$lib/utils/types";
import { error, redirect } from "@sveltejs/kit";
import type { PageServerLoad, RouteParams } from "./$types";


let api = new CerebroApi(private_env.PRIVATE_CEREBRO_API_URL_DOCKER);

const registerTrainingSession = async(fetch: Function, team_param: string, params: RouteParams): Promise<string | null> => {
    
    const registerSessionRequestInit: RequestInit = {
        method: 'POST',
        mode: 'cors',
        credentials: 'include',
        headers: { 'Content-Type': 'application/json' }, 
        body: JSON.stringify({
            "collection": params.collection,
            "shuffle": true
        } as CreateTrainingSession) 
    };

    let sessionIdentifier: string | null = null;

    let registerTrainingSessionResponse: Response = await fetch(
        `${api.routes.training.registerSession}?team=${team_param}`, registerSessionRequestInit
    );

    if (registerTrainingSessionResponse.ok) {
        
        // Store the training session identifier for the newly created session
        let registerTrainingSessionResponseData: TrainingResponse<string> = await registerTrainingSessionResponse.json();
        sessionIdentifier = registerTrainingSessionResponseData.data;

    } else {
        let errorResponse: ErrorResponseData = await registerTrainingSessionResponse.json();
        throw error(registerTrainingSessionResponse.status, errorResponse.message)
    }

    return sessionIdentifier
}

const fetchTrainingData = async(fetch: Function, team_param: string, params: RouteParams): Promise<[TrainingRecord, TrainingPrefetchData] | null> => {


    // Fetch the actual collection data for this session
    const fetchTrainingDataRequestInit: RequestInit = {
        method: 'GET',
        mode: 'cors',
        credentials: 'include'
    };

    let data: [TrainingRecord, TrainingPrefetchData] | null = null

    let fetchTrainingDataResponse: Response = await fetch(
        `${api.routes.training.getData}/${params.session}?team=${team_param}&record=${params.record}`, fetchTrainingDataRequestInit
    );

    if (fetchTrainingDataResponse.ok) {   
        let fetchTrainingDataResponseData: TrainingResponse<[TrainingRecord, TrainingPrefetchData]> = await fetchTrainingDataResponse.json();
        data = fetchTrainingDataResponseData.data;
    } else {
        let errorResponse: ErrorResponseData = await fetchTrainingDataResponse.json();
        throw error(fetchTrainingDataResponse.status, errorResponse.message)
    }

    return data
}

export const load: PageServerLoad = async ({ params, locals, fetch, depends, parent }) => { 

    depends("training:data")

    const { selectedTeam } = await parent();

    let data: [TrainingRecord, TrainingPrefetchData] | null = null;

    // Special case on initial load of the page for default values
    if (params.session === "0") {
        // Start a new session by registering it with the team database
        let sessionIdentifier = await registerTrainingSession(fetch, selectedTeam.name, params);
        throw redirect(302, `/cerebro/training/team=${selectedTeam.name}/collection=${params.collection}&session=${sessionIdentifier}&record=0`)
    } else {
        // Try and load training data from the team database
        data = await fetchTrainingData(fetch, selectedTeam.name, params)
    }

    return { 
        trainingData: data ? data[1] : null,
        trainingRecord: data ? data[0] : null
    };

}