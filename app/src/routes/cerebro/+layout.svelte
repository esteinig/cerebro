<script lang="ts">

import { Modal, ProgressBar, initializeStores } from '@skeletonlabs/skeleton';
import { AppShell } from '@skeletonlabs/skeleton';
import { Toast } from '@skeletonlabs/skeleton';

import Header from '$lib/components/shell/Header.svelte';
import Footer from '$lib/components/shell/Footer.svelte';
import { navigating } from '$app/stores';

initializeStores();

</script>

<Modal />
<Toast position="bl" color="variant-filled-tertiary"/>

<AppShell>
	<svelte:fragment slot="header">
        <Header />
		{#if $navigating}
			<ProgressBar height="h-1"/>
		{/if}
	</svelte:fragment>
    
	<svelte:fragment slot="footer">
        <Footer />
	</svelte:fragment>

    <slot />

</AppShell>