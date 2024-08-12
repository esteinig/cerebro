import type { Team, TeamDatabase } from "$lib/utils/types";
import { error, redirect } from "@sveltejs/kit";
import type { PageServerLoad } from "./$types";

export const load: PageServerLoad = async ({ locals, url, }) => {

    let currentUserTeams: Team[] = locals.teams;

    // Pick the first team and database from the user 
    // return an error if they are not part of a team
    // (which always has one default database)

    if (!currentUserTeams.length){
        throw error(404, `You are not a member of any teams - ${locals.admin ? "create a team for data upload first." : "please contact your administrator."}`)
    }

    let team: Team = currentUserTeams[0];

    throw redirect(302, `${url}/team=${team.id}&watcher=0`)

}