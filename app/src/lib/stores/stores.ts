import { browser } from "$app/environment";
import { writable, type Writable } from "svelte/store";

// Session based theme store. Grabs the current theme from the current body.
export const storeTheme = writable(browser ? document.body.getAttribute('data-theme') ?? '' : 'wintry');


interface SessionSettings {
    acceptedCommentRisk: boolean
}

// Session based settings store. Keeps user settings like risk acceptance checks.
export const storeSettings: Writable<SessionSettings> = writable({
    acceptedCommentRisk: false
} satisfies SessionSettings);

