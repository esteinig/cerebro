<div class="container mx-auto text-left p-4 md:text-center">
	<div class="relative z-[1] space-y-6 pt-32 pb-24">
		<div class="font-bold text-4xl tracking-tight">
			<span class="text-primary-500">Cerebro</span> <span class="text-3xl">metagenomics</span>
		</div>
		<p class="max-w-lg opacity-70 md:mx-auto">
			Diagnostic service for <span class="text-secondary-400">central nervous system</span> infections
			at the Victorian Infectious Disease Reference Laboratory (VIDRL)
		</p>
		<div class="flex flex-wrap gap-3 md:space-x-3 md:block">
			<a href="/login" class="btn md:btn-lg w-full md:w-fit variant-ringed-primary">Login</a>
			<a href="https://docs.meta-gp.org" class="btn md:btn-lg w-full md:w-fit variant-ringed-secondary">Docs</a>
		</div>
	</div>
</div>
<script>
	const mockdata = [

		{
			title: 'Clinical',
			description: 'Cerebro is used as a diagnostic service for life-threatening neurological infections in the state of Victoria, Australia.'
		},
		{
			title: 'Accreditation',
			description: 'Cerebro and the main sequencing assay are NATA-accredited at the Victorian Infectious Disease Reference Laboratory (VIDRL).'
		},
		{
			title: 'Development',
			description: "Cerebro is developed as part of a national metagenomics diagnostic project (Meta-GP) funded by the Medical Research Future Fund (MRFF)."
		},
		{
			title: 'Analysis',
			description: 'Assay protocol and analysis pipeline implement best practices in the field of metagenomic pathogen detection. Our work is published in ...'
		},
		{
			title: 'Open-source',
			description: 'Code and assay protocols are open-source under the GNU GPL v3. Cerebro can be deployed in local networks for increased data privacy and security.',
		},
		{
			title: 'Research',
			description: 'Cerebro is implemented by a multi-disciplinary clinical research team at The Peter Doherty Institute for Infection and Immunity with collaborations at the UK Health Security Agency and the Mayo Clinic.'
		},
	]
</script>

<div class="container mx-auto max-w-screen-lg">
	<div class="text-center space-y-4">
		<p class="max-w-xl opacity-60 mx-auto text-sm md:text-lg">
			Cerebro is a full-stack application written in <span class="text-secondary-400">Rust</span> and <span class="text-secondary-400">Svelte</span> driven by metagenomic pathogen identification in <span class="text-secondary-400">Nextflow</span>. 
		</p>
	</div>
	<div class="grid grid-cols-1 sm:grid-cols-2 md:grid-cols-3 gap-4 pt-16">
		{#each mockdata as feature}
			<div class="p-4 space-y-2 text-center">

				<div class="text-center">
					{feature.title}
				</div>
				<div class="text-xs opacity-60">
					{feature.description}
				</div>
			</div>
		{/each}
	</div>
	<div class="flex justify-center mt-10">
		<a href="https://github.com/esteinig/cerebro" type="button" class="btn-icon btn-icon-sm variant-filled">
			<svg aria-hidden="true" class="octicon octicon-mark-github" height="24" version="1.1" viewBox="0 0 16 16" width="24">
				<path fill-rule="evenodd" d="M8 0C3.58 0 0 3.58 0 8c0 3.54 2.29 6.53 5.47 7.59.4.07.55-.17.55-.38 0-.19-.01-.82-.01-1.49-2.01.37-2.53-.49-2.69-.94-.09-.23-.48-.94-.82-1.13-.28-.15-.68-.52-.01-.53.63-.01 1.08.58 1.23.82.72 1.21 1.87.87 2.33.66.07-.52.28-.87.51-1.07-1.78-.2-3.64-.89-3.64-3.95 0-.87.31-1.59.82-2.15-.08-.2-.36-1.02.08-2.12 0 0 .67-.21 2.2.82.64-.18 1.32-.27 2-.27.68 0 1.36.09 2 .27 1.53-1.04 2.2-.82 2.2-.82.44 1.1.16 1.92.08 2.12.51.56.82 1.27.82 2.15 0 3.07-1.87 3.75-3.65 3.95.29.25.54.73.54 1.48 0 1.07-.01 1.93-.01 2.2 0 .21.15.46.55.38A8.013 8.013 0 0 0 16 8c0-4.42-3.58-8-8-8z"></path>
			</svg>
		</a>
	</div>
</div>


