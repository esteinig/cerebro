import type { PageServerLoad } from './$types';
import { postJSON } from '$lib/server/api';
import { env as private_env } from "$env/dynamic/private";
import { redirect } from '@sveltejs/kit';

type VerificationCheckOK = { status: string; message: string; access_token: string };

export const load: PageServerLoad = async ({ url, fetch }) => {

  // If redirected here after success
  if (url.searchParams.get('complete') === 'true') {
    return { step: 'email_sent' as const };
  }

  const token = url.searchParams.get('token');

  // No token? Show instructions.
  if (!token) {
    return { step: 'needs_token' as const };
  }


  // Exchange email link token for a short-lived PasswordResetEmail token.
  // Immediately trigger the email send with that token.
  try {
    const check = await postJSON<VerificationCheckOK>('/auth/verification-check', {
      access_token: token
    }, private_env.PRIVATE_CEREBRO_API_URL_DOCKER);

    // Now request the password reset email be sent.
    await postJSON('/auth/password-reset-email', {
      access_token: check.access_token
      }, private_env.PRIVATE_CEREBRO_API_URL_DOCKER);

  } catch {
    // Treat any error as invalid/expired
    return { step: 'invalid' as const };
  }
  
  // Remove token, add complete=true, then redirect on successful password reset
  let clean = new URL(url);

  clean.searchParams.delete('token');
  clean.searchParams.set('complete', 'true');

  // Prevents form re-submits and avoids duplicate loads
  throw redirect(303, clean.pathname + (clean.searchParams.toString() ? `?${clean.searchParams.toString()}` : ''));

};