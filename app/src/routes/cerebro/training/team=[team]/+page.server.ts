import { env as private_env } from "$env/dynamic/private";
import CerebroApi from "$lib/utils/api";
import type { Team } from "$lib/utils/types";
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
    
    return { 
        selectedTeam: team
    };

}