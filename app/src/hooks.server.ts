
import { error, redirect } from '@sveltejs/kit';
import { env as public_env } from '$env/dynamic/public';
import { env as private_env } from '$env/dynamic/private';
import { parseSameSite, parseBool } from '$lib/server/api';

import {
    Role,
    type UserResponseData, 
    type AuthRefreshErrorResponseData,
    type AuthRefreshResponseData,
    type UserTeamsResponseData,
} from '$lib/utils/types';

import { CerebroApi } from '$lib/utils/api';

let api = new CerebroApi(private_env.PRIVATE_CEREBRO_API_URL_DOCKER);

// Handle all route requests and place route guards based on user role

/** @type {import('@sveltejs/kit').Handle} */
export async function handle({ event, resolve }) {
    
    const requestedPath = event.url.pathname;
    const cookies = event.cookies;

    let hookExceptions = requestedPath === "/" || 
        requestedPath.startsWith('/login') ||
        requestedPath.startsWith('/password') ||
        requestedPath.startsWith('/verify') ||    
        requestedPath.startsWith('/public')
        
    if (hookExceptions){
        // If the user requests the home, login, account confirmation, 
        // or password reset pages, or any public pages, don't modify 
        // the request and return response
        const response = await resolve(event, {
            transformPageChunk: ({ html }) => html.replace('%skeletonTheme%', cookies.get("skeletonTheme") ?? "dali"),
        });
        return response;
    }

    let accessToken: string | undefined = cookies.get("access_token");
    
    let userRequestInitOptions: RequestInit = {
        method: 'GET',
        mode: 'cors',
        credentials: 'include',
        headers: { cookie: `access_token=${accessToken}` }
    };    


    // this fetch is done only to check for access token validity
    // we again fetch user and their team data below in case
    // we need to refresh tokens and repeat check
    let userResponse: Response = await fetch(
        api.routes.users.self, userRequestInitOptions
    )

    if (!userResponse.ok) {
        // A failed response from the authentication middleware with potential indication
        // of possibility of refreshing access token (expired)
        let userResponseErrorData: AuthRefreshErrorResponseData = await userResponse.json();

        // Try and refresh the token, if token expired
        // explicitly set the refresh cookie here as we
        // are server-side and credentials of non same-site
        // origin are not passed on (see below)

        if (userResponseErrorData.refresh) {
            
            let refresh_token = cookies.get("refresh_token")

            const refreshRequestInitOptions: RequestInit = {
                method: 'GET',
                mode: 'cors',
                credentials: 'include',
                headers: { cookie: `refresh_token=${refresh_token}` }
            };

            let refreshResponse: Response = await fetch(
                api.routes.auth.refresh, 
                refreshRequestInitOptions
            )
            
            if (refreshResponse.ok){
                let refreshResponseData: AuthRefreshResponseData = await refreshResponse.json()

                // Explicitly reset the access token, since server-side
                cookies.set('access_token', refreshResponseData.access_token, {
                    path: '/',
                    httpOnly: true,
                    domain: private_env.PRIVATE_CEREBRO_API_ACCESS_COOKIE_DOMAIN,
                    secure: parseBool(private_env.PRIVATE_CEREBRO_API_ACCESS_COOKIE_SECURE),          // check this on insecure deployment (http that is not localhost)
                    sameSite: parseSameSite(private_env.PRIVATE_CEREBRO_API_ACCESS_COOKIE_SAME_SITE),                               // check this on insecure deployment (http that is not localhost)
                    maxAge: parseInt(private_env.PRIVATE_CEREBRO_API_ACCESS_MAX_AGE) * 60
                });

                // Pass on the now valid access token
                accessToken = refreshResponseData.access_token;

                userRequestInitOptions = {
                    method: 'GET',
                    mode: 'cors',
                    credentials: 'include',
                    headers: { cookie: `access_token=${accessToken}` }
                }

                userResponse = await fetch(
                    api.routes.users.self, userRequestInitOptions
                )

            } else {
                // If refresh request fails, redirect to login
                throw redirect(307, '/login')
            }
        } else {
            throw error(userResponse.status, userResponseErrorData.message)
        }
    }

    
    // Parallel requests to the user self routes to get user data and teams
    const userSelfTeamsResponse: Response= await fetch(api.routes.users.selfTeams, userRequestInitOptions)

    if (userResponse.ok && userSelfTeamsResponse.ok){
        let userSelfResponseData: UserResponseData = await userResponse.json();
        let userSelfTeamsResponseData: UserTeamsResponseData = await userSelfTeamsResponse.json();

        // Assign data to locals for availability on server-side load functions
        event.locals.admin = userSelfResponseData.data.user.roles.includes(Role.Admin);
        event.locals.user = userSelfResponseData.data.user;
        event.locals.teams = userSelfTeamsResponseData.data.teams;
    } else {
        throw error(500, "An unknown error during user data access occurred")
    }


    // Finally check if the user role allows for protected access to the routes
    if (requestedPath.includes("/admin") && !event.locals.user.roles.includes(Role.Admin)) {
        throw error(401, "Not authorized");
    } else if (requestedPath.includes("/cerebro") && !event.locals.user.roles.includes(Role.User)){
        throw error(401, "Not authorized");
    }

    let response = await resolve(event, {
        transformPageChunk: ({ html }) => html.replace('%skeletonTheme%', cookies.get("skeletonTheme") ?? "dali"),
    });

    return response;
}



/** @type {import('@sveltejs/kit').HandleFetch} */
export async function handleFetch({ event, request, fetch }) {

    if (request.url.startsWith(public_env.PUBLIC_CEREBRO_API_URL)) {
        // Clone the original request, but change the URL
        request = new Request(
            request.url.replace(public_env.PUBLIC_CEREBRO_API_URL, private_env.PRIVATE_CEREBRO_API_URL_DOCKER),  // change to docker container address on container network
            request
        );
    }

    
    request.headers.set('cookie', `refresh_token=${event.cookies.get("refresh_token")}; access_token=${event.cookies.get("access_token")}`);

    return fetch(request);
}
