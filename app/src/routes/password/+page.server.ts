import type { Actions, PageServerLoad } from './$types';
import { postJSON, parseBool, parseSameSite } from '$lib/server/api';
import { redirect, fail } from '@sveltejs/kit';
import { env as private_env } from "$env/dynamic/private";

type ResetCheck =
  | { status: 'ok'; access_token: string; message?: string }
  | { status: 'already_used'; access_token?: string; message?: string };

const COOKIE = 'reset_token';

export const load: PageServerLoad = async ({ url, cookies }) => {
  const token = url.searchParams.get('token');

  if (token) {
    
    let check: ResetCheck;
    try {
      check = await postJSON<ResetCheck>(
        '/auth/password-reset-check',
        { access_token: token },
        private_env.PRIVATE_CEREBRO_API_URL_DOCKER
      );
    } catch {
      return { step: 'invalid_request' as const };
    }

    if (check.status === 'ok' || (check.status === 'already_used' && check.access_token)) {
      cookies.set(COOKIE, check.access_token!, {
        httpOnly: true,
        sameSite: parseSameSite(private_env.PRIVATE_CEREBRO_API_ACCESS_COOKIE_SAME_SITE ?? 'strict'),
        secure: parseBool(private_env.PRIVATE_CEREBRO_API_ACCESS_COOKIE_SECURE ?? 'true'),
        domain: private_env.PRIVATE_CEREBRO_API_ACCESS_COOKIE_DOMAIN,
        path: '/password',
        maxAge: 10 * 60
      });
      throw redirect(303, '/password');
    }

    return { step: 'needs_token' as const };
  }

  const hasCookie = Boolean(cookies.get(COOKIE));
  return hasCookie ? { step: 'form' as const } : { step: 'needs_token' as const };
};

export const actions: Actions = {
  reset: async ({ request, cookies }) => {
    const form = await request.formData();
    const password = String(form.get('password') ?? '');
    const confirm = String(form.get('confirm') ?? '');

    if (password.length < 12 || password.length > 128) return fail(400, { error: 'Password length must be 12–128.' });
    if (password !== confirm) return fail(400, { error: 'Passwords do not match.' });

    const access_token = cookies.get(COOKIE);
    if (!access_token) return fail(403, { error: 'Reset session missing or expired.' });

    try {
      await postJSON('/auth/password-reset', { access_token, password }, private_env.PRIVATE_CEREBRO_API_URL_DOCKER);
    } finally {
      cookies.delete(COOKIE, { path: '/password' });
    }

    throw redirect(303, '/login');
  }
};