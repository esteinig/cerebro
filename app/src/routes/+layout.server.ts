/** @type {import('./$types').LayoutServerLoad} */
export async function load({ fetch, cookies, locals }) {
    return { 
        theme: cookies.get("skeletonTheme") ?? "wintry"
     }
}