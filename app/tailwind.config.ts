import { join } from 'path';
import type { Config } from 'tailwindcss';
import forms from '@tailwindcss/forms';
import typography from '@tailwindcss/typography';
import { skeleton } from '@skeletonlabs/tw-plugin';
import { momaDali } from './themes/dali';

export default {
	darkMode: 'class',
	content: ['./src/**/*.{html,js,svelte,ts}', join(require.resolve('@skeletonlabs/skeleton'), '../**/*.{html,js,svelte,ts}')],
	theme: {
		extend: {},
	},
	plugins: [
		forms,
		typography,
		skeleton({
			themes: {
				custom: [
					momaDali
				],
				preset: [
					{
						name: 'hamlindigo',
						enhancements: true,
					},
					{
						name: 'wintry',
						enhancements: true,
					},
					{
						name: 'crimson',
						enhancements: true,
					},
					{
						name: 'skeleton',
						enhancements: true,
					},
					{
						name: 'modern',
						enhancements: true,
					},
					{
						name: 'rocket',
						enhancements: true,
					},
					{
						name: 'seafoam',
						enhancements: true,
					},
					{
						name: 'vintage',
						enhancements: true,
					},
					{
						name: 'sahara',
						enhancements: true,
					},
					{
						name: 'gold-nouveau',
						enhancements: true,
					}
				],
			},
		}),
	],
} satisfies Config;