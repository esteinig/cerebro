<script lang="ts">
    import { page } from "$app/stores";
</script>



<div class="card max-w-xl2">
	<header class="card-header">
        <div class="p-4">
            <div class="flex">
                <div class="flex items-center">
                    <div class="flex-none ml-5">
                        <div class="md:h4 lg:h3">{$page.data.userData.title ?? ""} {$page.data.userData.name}</div>
                        <div class="flex items-center">
                            <h4 class="h4 opacity-60 pr-1">{$page.data.userData.positions.join(", ")}</h4>
                        </div>
                    </div>
                    

                </div>
            </div>
        </div>
    </header>
    <section>

    </section>
	<footer class="card-footer">

        
    </footer>
</div>