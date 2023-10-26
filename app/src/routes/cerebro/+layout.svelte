<script lang="ts">

import { Modal, initializeStores } from '@skeletonlabs/skeleton';
import { AppShell } from '@skeletonlabs/skeleton';
import { Toast } from '@skeletonlabs/skeleton';

import Header from '$lib/components/shell/Header.svelte';
import Footer from '$lib/components/shell/Footer.svelte';

initializeStores();

</script>

<Modal />
<Toast position="bl" color="variant-filled-tertiary"/>

<AppShell>
	<svelte:fragment slot="header">
        <Header />
	</svelte:fragment>
    
	<svelte:fragment slot="footer">
        <Footer />
	</svelte:fragment>

    <slot />

</AppShell>