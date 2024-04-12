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

    if (!team.databases.length){
        throw error(404, `No database has been created for this team (${team.name})`)
    }

    let teamDatabase: TeamDatabase = team.databases[0];

    throw redirect(302, `${url}/team=${team.id}&db=${teamDatabase.id}&location=0`)

}