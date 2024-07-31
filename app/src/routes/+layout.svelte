<script lang="ts">

	import '../app.postcss';
	
	// Highlight JS
	import hljs from 'highlight.js';
	import 'highlight.js/styles/github-dark.css';
	import { storeHighlightJs } from '@skeletonlabs/skeleton';
	storeHighlightJs.set(hljs);
	
	// Floating UI for Popups
	import { computePosition, autoUpdate, flip, shift, offset, arrow } from '@floating-ui/dom';
	import { storePopup } from '@skeletonlabs/skeleton';
	import { browser } from '$app/environment';
	import { storeTheme } from '$lib/stores/stores';

	export let data;

	storePopup.set({ computePosition, autoUpdate, flip, shift, offset, arrow });

	$storeTheme = data.theme ?? "%skeletonTheme%";

	const setBodyThemeAttribute = (): void => {
		if (!browser) return;
		document.body.setAttribute('data-theme', $storeTheme);
	}

	storeTheme.subscribe(setBodyThemeAttribute);
	

</script>


<slot />
