import type { PageServerLoad } from './$types';
import { postJSON } from '$lib/server/api';
import { env as private_env } from "$env/dynamic/private";
import { redirect } from '@sveltejs/kit';


type VerificationCheckOK = { status: string; message: string; access_token: string };

export const ssr = true;
export const csr = false;

const COOKIE = 'vtoken';

// in-process dedupe across concurrent requests in this instance
const inflight = new Map<string, Promise<void>>();

async function verifyAndSend(token: string) {
  const check = await postJSON<VerificationCheckOK>(
    '/auth/verification-check',
    { access_token: token },
    private_env.PRIVATE_CEREBRO_API_URL_DOCKER
  );
  await postJSON(
    '/auth/password-reset-email',
    { access_token: check.access_token },
    private_env.PRIVATE_CEREBRO_API_URL_DOCKER
  );
}

export const load: PageServerLoad = async ({ url, cookies }) => {
  const token = url.searchParams.get('token');
  const step = url.searchParams.get('step');

  if (step === 'done') return { step: 'email_sent' as const };

  // First hop: stash token, move to clean URL
  if (token) {
    cookies.set(COOKIE, token, { path: '/', httpOnly: true, sameSite: 'lax', secure: true, maxAge: 300 });
    const u = new URL(url);
    u.searchParams.delete('token');
    u.searchParams.set('step', 'process');
    throw redirect(303, `${u.pathname}?${u.searchParams}`);
  }

  if (step === 'process') {
    const handed = cookies.get(COOKIE);
    if (!handed) return { step: 'invalid' as const };

    // dedupe: share one in-flight call per token
    let p = inflight.get(handed);
    if (!p) {
      p = verifyAndSend(handed).finally(() => inflight.delete(handed));
      inflight.set(handed, p);
    }

    try {
      await p;
    } catch {
      cookies.delete(COOKIE, { path: '/' });
      return { step: 'invalid' as const };
    }

    cookies.delete(COOKIE, { path: '/' });
    const u = new URL(url);
    u.searchParams.set('step', 'done');
    throw redirect(303, `${u.pathname}?${u.searchParams}`);
  }

  return { step: 'needs_token' as const };
};