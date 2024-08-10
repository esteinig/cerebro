<script lang="ts">
    export let active: boolean = false;
    export let size: string = "4";
    export let divClass: string = "";
</script>

<div
    class={`${divClass} rounded-full ${active ? 'bg-green-500 animate-blink' : 'bg-red-500'}`}
    style={`width: ${size}rem; height: ${size}rem;`}
></div>

<style>
    @keyframes blink {
        0%, 100% {
            opacity: 1;
        }
        50% {
            opacity: 0;
        }
    }

    .animate-blink {
        animation: blink 3s infinite;
    }
</style>