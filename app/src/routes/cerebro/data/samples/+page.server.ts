import { env as private_env } from "$env/dynamic/private";
import CerebroApi from "$lib/utils/api";
import type { ProjectCollection, SampleOverviewData, SampleOverviewResponseData, Team, TeamDatabase } from "$lib/utils/types";
import { error } from "@sveltejs/kit";
import type { PageServerLoad } from "./$types";

let api = new CerebroApi(private_env.PRIVATE_CEREBRO_API_URL_DOCKER);

export const load: PageServerLoad = async ({ locals, fetch }) => {

    let currentUserTeams: Array<Team> = locals.teams;

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

    // Pick the first project from the team database 
    // return an error if no projects have been 
    // registered for this team database

    if (!teamDatabase.projects.length){
        throw error(404, `No projects have been registered for this database (${teamDatabase.name})`);
    }

    let project: ProjectCollection = teamDatabase.projects[0];

    // Here we can get the negative template control tag 
    // and page limits from the user settings later

    const negativeTemplateControl: string = "NTC";
    const pageLimit: number = 500;

    // Fetch the sample overview data
    const sampleOverviewResponse: Response = await fetch(
        `${api.routes.cerebro.sampleOverview}?team=${team.id}&db=${teamDatabase.id}&project=${project.id}&page=0&limit=${pageLimit}`, 
        { method: 'GET', mode: 'cors', credentials: 'include' } as RequestInit
    );

    // Can be empty if none found so that another team or project can still be selected
    let sampleOverviewData: Array<SampleOverviewData> = [];
    
    if (sampleOverviewResponse.ok) {
        let sampleOverviewResponseData: SampleOverviewResponseData = await sampleOverviewResponse.json();

        // Do not thow an error when no samples have been found in this database / project
        // so that the page loads - this is a valid response

        sampleOverviewData = sampleOverviewResponseData.data.sample_overview;

    } else {
        let responseErrorData = await sampleOverviewResponse.json();
        throw error(sampleOverviewResponse.status, responseErrorData.message)
    }
    
    return { 
        sampleOverview: sampleOverviewData,
        defaultTeam: team,
        defaultDatabase: teamDatabase,  
        defaultProject: project, 
        defaultNegativeTemplateControl: negativeTemplateControl, 
        defaultPageLimit: pageLimit 
    };

    
}