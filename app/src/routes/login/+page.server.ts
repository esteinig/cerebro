import type { Actions } from './$types';

export const actions: Actions = {
	// This action is called when the user clicks the login button
	login: async ({ cookies, request }) => {

		const formData = await request.formData();
		const theme = formData.get('theme')?.toString() ?? 'dali';
		
		// Sets the selected theme to the cookie
		cookies.set('theme', theme, { path: '/' });
		return { theme };
	}
};