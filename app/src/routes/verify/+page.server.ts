import type { PageServerLoad } from './$types';
import { postJSON } from '$lib/server/api';
import { env as private_env } from "$env/dynamic/private";
import { redirect } from '@sveltejs/kit';

type VerificationCheckOK = { status: string; message: string; access_token: string };

export const ssr = true;
export const csr = false;

export const load: PageServerLoad = async ({ url, cookies }) => {
  const token = url.searchParams.get('token');
  const step = url.searchParams.get('step');

  // step=done renders success without side effects
  if (step === 'done') return { step: 'email_sent' as const };

  // First hit with token: stash to cookie and move to a clean URL
  if (token) {
    cookies.set('vtoken', token, { path: '/', httpOnly: true, sameSite: 'lax', secure: true, maxAge: 300 });
    const clean = new URL(url);
    clean.searchParams.delete('token');
    clean.searchParams.set('step', 'process');
    throw redirect(303, `${clean.pathname}?${clean.searchParams}`);
  }

  // Only this URL performs side effects. Bots that just GET the token URL won’t reach here.
  if (step === 'process') {
    
    const handed = cookies.get('vtoken');
    if (!handed) return { step: 'invalid' as const };

    try {
      const check = await postJSON<VerificationCheckOK>(
        '/auth/verification-check',
        { access_token: handed },
        private_env.PRIVATE_CEREBRO_API_URL_DOCKER
      );

      await postJSON(
        '/auth/password-reset-email',
        { access_token: check.access_token },
        private_env.PRIVATE_CEREBRO_API_URL_DOCKER
      );

      cookies.delete('vtoken', { path: '/' });
      const done = new URL(url);
      done.searchParams.set('step', 'done');
      throw redirect(303, `${done.pathname}?${done.searchParams}`);
    } catch {
      cookies.delete('vtoken', { path: '/' });
      return { step: 'invalid' as const };
    }
  }

  // No token and not processing
  return { step: 'needs_token' as const };
};