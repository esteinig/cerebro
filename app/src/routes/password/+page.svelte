<script lang="ts">
    export let data:
      | { step: 'needs_token' }
      | { step: 'request_new' }
      | { step: 'form' }
      | { step: 'invalid' };
    export let form:
      | { error?: string }
      | undefined;
    let show = false;
</script>

<div class="container h-full mx-auto flex justify-center items-center">

    {#if data.step === 'needs_token'}
        <p class="h2">Open the reset link from your email</p>
    {/if}
    {#if data.step === 'request_new'}
        <p class="h2">Please request a new link</p>
    {/if}
    {#if data.step === 'invalid'}
        <p class="h2">Link expired or invalid</p>
    {/if}

    {#if data.step === 'form'}
        <div class="card p-6 space-y-6 shadow-xl text-left"  data-sveltekit-preload-data="off">
            
            <form method="POST" class="space-y-4">
                <div>
                    <label class="label" for="password">New password</label>
                    <input class="input" id="password" name="password" type={show ? 'text' : 'password'} minlength="12" maxlength="128" required autocomplete="new-password" />
                </div>
                <div>
                    <label class="label" for="confirm">Confirm password</label>
                    <input class="input" id="confirm" name="confirm" type={show ? 'text' : 'password'} minlength="12" maxlength="128" required autocomplete="new-password" />
                </div>
                <label class="label" >
                    <input class="checkbox" type="checkbox" bind:checked={show} />
                    Show passwords
                </label>
                <button class="btn variant-filled-primary w-full" formaction="?/reset" formmethod="POST">Reset password</button>
            </form>
            {#if form?.error}
                <div class="space-y-2">
                    <div class="flex justify-center">
                        <svg aria-hidden="true" fill="none" stroke="currentColor" class="h-16 w-16 text-tertiary-500" stroke-width="1" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg">
                            <path d="M12 9v3.75m9-.75a9 9 0 11-18 0 9 9 0 0118 0zm-9 3.75h.008v.008H12v-.008z" stroke-linecap="round" stroke-linejoin="round"></path>
                        </svg> 
                    </div>
                    <p class="text-center" aria-live="polite">{form.error}</p>  
                </div> 
            {/if}
        </div>
    {/if}
</div>
  