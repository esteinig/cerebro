
<script lang="ts">
	import CheckmarkIcon from "../icons/CheckmarkIcon.svelte";
	import CircleIndicator from "../icons/CircleIndicator.svelte";

    export let chips: string[];
    export let divClass: string = "";
    export let circle: boolean = false;
    export let circleColor: string = "bg-primary-500";
    export let activeCondition: boolean = false;
    export let clickFunction: (chip: string, condition: boolean) => void = () => {};

</script>

<div class={divClass}>
    {#each chips as chip}
        <span class="chip {activeCondition ? 'variant-ghost-primary' : 'variant-soft'} mr-2" on:click={() => clickFunction?.(chip, !activeCondition)} on:keypress aria-hidden>
            {#if activeCondition}
                <span><CheckmarkIcon/></span>
            {/if}
            <span>{chip}</span>
            {#if circle}
                <CircleIndicator color={circleColor}/>
            {/if}
        </span>
    {/each}
</div>
   