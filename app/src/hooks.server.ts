
import { error, redirect } from '@sveltejs/kit';
import { env as public_env } from '$env/dynamic/public';
import { env as private_env } from '$env/dynamic/private';

import {
    Role,
    type UserResponseData, 
    type AuthRefreshErrorResponseData,
    type AuthRefreshResponseData,
    type UserTeamsResponseData,
    type ErrorResponseData 
} from '$lib/utils/types';

import { CerebroApi } from '$lib/utils/api';

let api = new CerebroApi(private_env.PRIVATE_CEREBRO_API_URL_DOCKER);

// Handle all route requests and place route guards based on user role

/** @type {import('@sveltejs/kit').Handle} */
export async function handle({ event, resolve }) {
    
    const requestedPath = event.url.pathname;
    const cookies = event.cookies;

    let hookExceptions = requestedPath === "/" || 
        requestedPath === "/login" || 
        requestedPath === "/password" || 
        requestedPath === "/verify"
        
    if (hookExceptions){
        // If the user requests the home, login, account confirmation, 
        // or password reset pages, don't modify the request and 
        // return response
        const response = await resolve(event, {
            transformPageChunk: ({ html }) => html.replace('%skeletonTheme%', cookies.get("skeletonTheme") ?? "dali"),
        });
        return response;
    }

    // For all other routes get the access token
    // and verify the user, returning the user
    // data from the API
    let accessToken: string = cookies.get("access_token");
    
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
                    secure: private_env.PRIVATE_CEREBRO_API_ACCESS_COOKIE_SECURE === "true" ? true : false,          // check this on insecure deployment (http that is not localhost)
                    sameSite: private_env.PRIVATE_CEREBRO_API_ACCESS_COOKIE_SAME_SITE,                               // check this on insecure deployment (http that is not localhost)
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

    // We also need to make sure we still have the user data
    // including the admin routes in locals, so we need another
    // call to the user self endpoint after refreshing, sending 
    // the new access token and checking again 
    
    
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

    // All this is so we can use the validated current user data to check permissions on the route

    // Finally check if the user role allows for protected access to the routes
    if (requestedPath.includes("/admin") && !event.locals.user.roles.includes(Role.Admin)) {
        throw error(401, "Not authorized");
    } else if (requestedPath.includes("/cerebro") && !event.locals.user.roles.includes(Role.User)){
        throw error(401, "Not authorized");
    }

    // BUG in routing: if we want to specify a page with a search param in the router e.g.
    // 
    // /meta-gp/samples/[sample][x+3f]db=[db] e.g.
    // /meta-gp/samples/Test?db=test
    //
    // then the `event.url` splits the off the search string from the path correctly
    //
    // `event.url.pathname` = /meta-gp/samples/Test
    // `event.url.search` = ?db=test
    //
    // However the resolve of the event will use the `event.url.pathname` not navigate to
    // the page without the search params 
    //
    // Error: Not found: /meta-gp/samples/Test
    //
    // If we assign manually the correct `event.url.pathname` before the resolve, the 
    // search param indicator `?` is not resolved correctly but translated into 
    // URL safe format:
    //
    // event.url.pathname = 'event.url.pathname' + 'event.url.search'  --> /meta-gp/samples/Test%3Fdb=test
    //
    // If we try to replace the file route format[x+3f], x+3f (or even [%3F]) the term is not resolved e.g.
    //
    // if (event.url.search !== ''){
    //    event.url.pathname = event.url.pathname + '[x+3f]' + event.url.search.substring(1)
    // }

    let response = await resolve(event, {
        transformPageChunk: ({ html }) => html.replace('%skeletonTheme%', cookies.get("skeletonTheme") ?? "dali"),
    });

    return response;
}


// Server-side fetches replaced with the Docker network address - this is needed
// because in hybrid files (e.g. `page.ts` or `layout.ts`) requests against the
// data and authentication API are executed both client- and server-side so that
// each request requires its appropriate API base URL (i.e. against public or 
// Docker in-network address)

// Handle all server-side fetch calls (e.g. in `+page.server.ts` or `+layout.server.ts`)

/** @type {import('@sveltejs/kit').HandleFetch} */
export async function handleFetch({ event, request, fetch }) {

    if (request.url.startsWith(public_env.PUBLIC_CEREBRO_API_URL)) {
        // Clone the original request, but change the URL
        request = new Request(
            request.url.replace(public_env.PUBLIC_CEREBRO_API_URL, private_env.PRIVATE_CEREBRO_API_URL_DOCKER),  // change to docker container address on container network
            request
        );
    }

    // Credentials are not included automatically:

    // For same-origin requests, SvelteKit's fetch implementation will forward cookie and authorization headers unless the credentials option is set to "omit".
    // For cross-origin requests, cookie will be included if the request URL belongs to a subdomain of the app — for example if your app is on my-domain.com, 
    // and your API is on api.my-domain.com, cookies will be included in the request. If your app and your API are on sibling subdomains — www.my-domain.com 
    // and api.my-domain.com for example — then a cookie belonging to a common parent domain like my-domain.com will not be included, because SvelteKit has
    // no way to know which domain the cookie belongs to. In these cases you will need to manually include the cookie using handleFetch.

    // Since we currently have the BASE_URL_API and BASE_URL_API_DOCKER on different domains
    // cookies are not automatically included, so we set them again manually:

    // OMG the request header, which should after token refresh contain the renewed `access_token` as described 
    // in the documentation to be accessed by: `event.request.headers.get('cookie')` DOES NOT contain the 
    // refreshed `access_token`. When not inserting token refresh or after reload of error page after toke nrefresh
    // it DOES correctly contain the refreshed `access_token`.

    // Solution: we have to get the cookie from the cookies object, likely because that's where it is set in the `handle` refresh
    // logic (above) and we make a server-side fetch request immediately after the refresh, so the access_token is not provided
    // by the browser (forwarded in `+page.server.ts`) but instead server-side by the cookies object itself (if that makes sense)
    
    request.headers.set('cookie', `refresh_token=${event.cookies.get("refresh_token")}; access_token=${event.cookies.get("access_token")}`);

    return fetch(request);
}
