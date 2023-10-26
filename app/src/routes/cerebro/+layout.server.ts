import { Role } from "$lib/utils/types";

/** @type {import('./$types').LayoutServerLoad} */
export async function load({ fetch, cookies, locals, depends }) {
    
    depends("cerebro:data")

    return { 
        admin: locals.user.roles.includes(Role.Admin),
        accessToken: cookies.get("access_token"),
        refreshToken: cookies.get("refresh_token"),
        userData: locals.user,
        userTeams: locals.teams
    }
}