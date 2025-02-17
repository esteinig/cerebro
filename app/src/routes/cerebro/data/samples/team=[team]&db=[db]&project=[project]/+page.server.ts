import { env as private_env } from "$env/dynamic/private";
import CerebroApi from "$lib/utils/api";
import type { ProjectCollection, SampleOverviewData, SampleOverviewResponseData, Team, TeamDatabase } from "$lib/utils/types";
import { error } from "@sveltejs/kit";
import type { PageServerLoad } from "./$types";

let api = new CerebroApi(private_env.PRIVATE_CEREBRO_API_URL_DOCKER);

export const load: PageServerLoad = async ({ params, locals, fetch }) => { 

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

    let database: TeamDatabase;

    if (!team.databases.length){
        throw error(404, `No database has been created for this team (${team.name})`)
    }

    // Special case on initial load of the page for default values
    if (params.db === "0") {
        database = team.databases[0];
    } else {
        let selectedDatabase = team.databases.find(db => db.id === params.db || db.name === params.db);

        if (selectedDatabase) { 
            database = selectedDatabase 
        } else { 
            throw error(404, "Database not found") 
        }
    }

    // Pick the first project from the team database 
    // return an error if no projects have been 
    // registered for this team database

    if (!database.projects.length){
        throw error(404, `No projects have been registered for this database (${database.name})`);
    }

    let project: ProjectCollection;


    // Special case on initial load of the page for default values
    if (params.project === "0") {
        project = database.projects[0];
    } else {
        
        let selectedProject = database.projects.find(project => project.id === params.project || project.name === params.project);

        if (selectedProject) { 
            project = selectedProject 
        } else { 
            throw error(404, "Project not found") 
        }
    }

    // Here we can get the negative template control tag 
    // and page limits from the user settings later

    const negativeTemplateControl: string = "NTC";
    const pageLimit: number = 500;

    // Fetch the sample overview data
    const sampleOverviewResponse: Response = await fetch(
        `${api.routes.cerebro.sampleOverview}?team=${team.id}&db=${database.id}&project=${project.id}&page=0&limit=${pageLimit}`, 
        { method: 'GET', mode: 'cors', credentials: 'include' } as RequestInit
    );

    // Can be empty if none found so that another team or project can still be selected
    let sampleOverview: Array<SampleOverviewData> = [];
    
    if (sampleOverviewResponse.ok) {
        let sampleOverviewResponseData: SampleOverviewResponseData = await sampleOverviewResponse.json();

        // Do not thow an error when no samples have been found in this database / project
        // so that the page loads - this is a valid response!

        sampleOverview = sampleOverviewResponseData.data.sample_overview;

    } else {
        let responseErrorData = await sampleOverviewResponse.json();
        throw error(sampleOverviewResponse.status, responseErrorData.message)
    }
    
    return { 
        sampleOverview: sampleOverview,
        selectedTeam: team,
        selectedDatabase: database,  
        selectedProject: project, 
        negativeTemplatecontrol: negativeTemplateControl, 
        pageLimit: pageLimit 
    };

    
}