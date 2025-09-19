import type { PageServerLoad } from './$types';
import { postJSON } from '$lib/server/api';
import { env as private_env } from "$env/dynamic/private";
import { redirect } from '@sveltejs/kit';


type VerificationCheckOK = { status: string; message: string; access_token: string };

export const load: PageServerLoad = async ({ url, fetch }) => {

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

    return { step: 'email_sent' as const };

  } catch {
    // Treat any error as invalid/expired
    return { step: 'invalid' as const };
  }

};