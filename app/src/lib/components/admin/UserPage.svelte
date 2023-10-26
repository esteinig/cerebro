<script lang="ts">
    import UserTable from "./tables/UserTable.svelte";
    import NewUserCard from "./cards/NewUserCard.svelte";
	import { page } from "$app/stores";
	import type { User } from "$lib/utils/types";
	import { getToastStore, type ToastSettings } from "@skeletonlabs/skeleton";
	import UserCard from "./cards/UserCard.svelte";

    const toastStore = getToastStore();

    let newUserCard: boolean = true;
    let selectedRow: string[] = [];

    let selectedUser: User | null = null;

    const getSelectedUser = (selectedRow: string[]): User | null => {

        if (selectedRow.length > 0){
            let userId: string = selectedRow[0];

            let match = $page.data.users.filter(
                (user: User) => user.id === userId
            )[0]
            if (match === undefined) {
                const errorToast: ToastSettings = {
                    message: "Failed to select user",
                    background: 'variant-filled-warning',
                };
                toastStore.trigger(errorToast);
                return null
            }
            return match
        }
        return null
    }   
    $: {
        selectedUser = getSelectedUser(selectedRow);
        if (selectedUser !== null) {
            newUserCard = false;
        }
    }



</script>

<div class="grid grid-cols-1 sm:grid-cols-1 md:grid-cols-2 gap-16 pt-10">
    <div class="">
        <UserTable bind:selectedRow/>
        <div class="flex justify-center gap-4 mt-2">
            <button type="button" class="btn variant-outline-success" on:click={() => newUserCard = true}>New User</button>
        </div>
    </div>
    <div class="">
        {#if newUserCard}
            <NewUserCard />
        {:else}
            {#if selectedUser !== null}
                <UserCard bind:user={selectedUser} bind:selectedRow/>
            {/if}
        {/if}
    </div>
</div>