// See https://kit.svelte.dev/docs/types#app
// for information about these interfaces

import type CerebroApi from "$lib/utils/api";
import type { RequestLog, Team, User } from "$lib/utils/types";

declare namespace App {
	interface PageData {
		// /cerebro
		admin?: boolean,
		publicApi?: CerebroApi,
		userData?: User,
		userTeams?: Array<Team>,
		// /admin
		users?: Array<User>,
		teams?: Array<Teams>,
		logs?: Array<RequestLog>,
		criticalLogs?: Array<RequestLog>
		
	}
	// interface Error {}
	interface Locals {
		// Available thorugh `hook.server.ts`
		// for other server-side page load
		// access
		admin?: boolean,
		user?: User,
		teams?: Team[]
	}
	// interface Platform {}
}
