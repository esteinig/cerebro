import type { Actions, PageServerLoad } from './$types';
import { postJSON, parseBool, parseSameSite } from '$lib/server/api';
import { redirect, fail } from '@sveltejs/kit';
import { env as private_env } from "$env/dynamic/private";

type ResetCheckOK = { status: string; message: string; access_token: string };

const COOKIE = 'reset_token';

function cookieSecure() {
  return parseBool(private_env.PRIVATE_CEREBRO_API_ACCESS_COOKIE_SECURE ?? 'true');
}

function cookieSameSite() {
  return parseSameSite(private_env.PRIVATE_CEREBRO_API_ACCESS_COOKIE_SAME_SITE ?? 'strict')
}

export const load: PageServerLoad = async ({ url, cookies }) => {
  const token = url.searchParams.get('token');

  // If this is the landing via email link, exchange and store short-lived token in an httpOnly cookie.
  if (token) {

    let access_token: string;

    try {
      ({ access_token } = await postJSON<ResetCheckOK>('/auth/password-reset-check', {
        access_token: token
      }, private_env.PRIVATE_CEREBRO_API_URL_DOCKER));

    } catch {
      return { step: 'invalid' as const };
    }

    // HttpOnly so JS cannot read it. Tight scope and lifetime.
    cookies.set(COOKIE, access_token, {
        httpOnly: true,
        sameSite: cookieSameSite(),
        secure: cookieSecure(),
        domain: private_env.PRIVATE_CEREBRO_API_ACCESS_COOKIE_DOMAIN,
        path: '/password',
        maxAge: 10*60
    });

    // Remove token from the URL.
    throw redirect(303, '/password');
  }

  // If already have the cookie, show the form. Otherwise, user hasn't come from a valid link.
  const hasCookie = Boolean(cookies.get(COOKIE));
  return hasCookie ? { step: 'form' as const } : { step: 'needs_token' as const };
};

export const actions: Actions = {
  reset: async ({ request, cookies }) => {
    const form = await request.formData();
    const password = String(form.get('password') ?? '');
    const confirm = String(form.get('confirm') ?? '');

    if (password.length < 12 || password.length > 128) {
      return fail(400, { error: 'Password length must be 12–128.' });
    }
    if (password !== confirm) {
      return fail(400, { error: 'Passwords do not match.' });
    }

    const access_token = cookies.get(COOKIE);
    if (!access_token) {
      return fail(403, { error: 'Reset session missing or expired.' });
    }

    try {
      await postJSON('/auth/password-reset', { access_token, password }, private_env.PRIVATE_CEREBRO_API_URL_DOCKER);
      cookies.delete(COOKIE, { path: '/password' });
    } catch {
      cookies.delete(COOKIE, { path: '/password' });
      return fail(403, { error: 'Link expired or invalid. Request a new reset email.' });
    }

    throw redirect(303, '/login');
  }
};