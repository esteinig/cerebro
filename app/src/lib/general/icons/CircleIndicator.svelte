<script lang="ts">
    export let size: number = 2;
    export let color: string = "bg-primary-500"
    export let circleClass: string = "";
</script>

<div class="rounded-full {color} h-{size} w-{size} {circleClass}"></div>