<script lang="ts">

import { Modal, ProgressBar, initializeStores } from '@skeletonlabs/skeleton';
import { AppShell } from '@skeletonlabs/skeleton';
import { Toast } from '@skeletonlabs/skeleton';

import Header from '$lib/components/shell/Header.svelte';
import Footer from '$lib/components/shell/Footer.svelte';
import { navigating } from '$app/stores';
import { navigationLoading } from '$lib/stores/stores';

initializeStores();

</script>

<Modal />
<Toast position="bl" color="variant-filled-tertiary"/>

<AppShell>
	<svelte:fragment slot="header">
        <Header />
		<div class="h-1 { $navigating || $navigationLoading ? 'opacity-100' : 'opacity-0' } transition-opacity ease-in duration-300 delay-100">
			<ProgressBar height="h-full"/>
		</div>
	</svelte:fragment>
    
	<svelte:fragment slot="footer">
        <Footer />
	</svelte:fragment>

    <slot />

</AppShell>