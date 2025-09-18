<script lang="ts">
    export let data:
      | { step: 'needs_token' }
      | { step: 'email_sent' }
      | { step: 'invalid' };
</script>

<div class="mt-[10%] flex justify-center text-lg">
{#if data.step === 'needs_token'}
    <p class="h2">Open the verification link from your email</p>
{/if}

{#if data.step === 'email_sent'}
    <p class="h2">Verification successful - check your inbox for a password reset email</p>
{/if}

{#if data.step === 'invalid'}
    <p class="h2">Link expired or invalid</p>
{/if}
</div>
  