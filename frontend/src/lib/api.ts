import axios, { type AxiosRequestConfig } from 'axios'
import type {
  RunTrainingForm,
  UserForm,
  UserLogin,
  CodonTableForm,
  UserPasswordForm,
  UserProfileForm,
  UserRoleForm,
  SequenceComparatorForm,
  TunedSequence,
} from './interfaces'
import { store } from './store'

interface ApiError {
  code: number | null
  detail: string
}

type ApiResponse<T> = [T, null] | [null, ApiError]

const client = axios.create({
  baseURL: '/api',
  withCredentials: true, // Enable sending cookies with requests
})

// Interceptor to unify error handling
client.interceptors.response.use(
  response => response,
  error => {
    if (axios.isAxiosError(error) && error.response) {
      return Promise.reject({
        code: error.response.status,
        detail: error.response.data.detail,
      })
    } else {
      const message = 'Network error or request failed'
      console.error(message)
      return Promise.reject({
        code: null,
        detail: message,
      })
    }
  },
)

let isRefreshing = false
let refreshPromise: Promise<boolean> | null = null

async function refreshAccessToken(): Promise<boolean> {
  try {
    // Call the refresh endpoint to get new tokens
    await client.post('/auth/refresh')
    return true
  } catch (error) {
    console.error('Failed to refresh token:', error)
    return false
  }
}

async function makeRequest(config: AxiosRequestConfig, isRetry = false) {
  try {
    const response = await client.request(config)
    return [response.data, null]
  } catch (error) {
    const apiError = error as ApiError

    // If we get a 401 and haven't already retried, attempt to refresh the token
    if (apiError.code === 401 && !isRetry) {
      // If a refresh is already in progress, wait for it
      if (isRefreshing && refreshPromise) {
        const refreshSucceeded = await refreshPromise
        if (refreshSucceeded) {
          // Retry the original request after successful refresh
          return makeRequest(config, true)
        } else {
          // Refresh failed, clear user state
          store.emptyCurrentUser()
        }
      } else {
        // Start a new refresh
        isRefreshing = true
        refreshPromise = refreshAccessToken()

        const refreshSucceeded = await refreshPromise

        isRefreshing = false
        refreshPromise = null

        if (refreshSucceeded) {
          // Retry the original request with the new access token
          return makeRequest(config, true)
        } else {
          // Refresh failed, clear user state
          store.emptyCurrentUser()
        }
      }
    }

    return [null, error] as [null, ApiError]
  }
}

async function login(form: UserLogin) {
  // Submit traditional form to allow browser to follow redirect and set cookie
  const formElement = document.createElement('form')
  formElement.method = 'POST'
  formElement.action = '/api/auth/token'

  const usernameInput = document.createElement('input')
  usernameInput.type = 'hidden'
  usernameInput.name = 'username'
  usernameInput.value = form.username

  const passwordInput = document.createElement('input')
  passwordInput.type = 'hidden'
  passwordInput.name = 'password'
  passwordInput.value = form.password

  formElement.appendChild(usernameInput)
  formElement.appendChild(passwordInput)

  document.body.appendChild(formElement)
  formElement.submit()

  // Form submission causes page reload, so this never returns
  return null
}

const REQUESTS = {
  get: async (url: string, params?: object) =>
    await makeRequest({ method: 'get', url, params }),
  post: async (url: string, data?: object) =>
    await makeRequest({ method: 'post', url, data }),
  put: async (url: string, data?: object) =>
    await makeRequest({ method: 'put', url, data }),
  delete: async (url: string) => await makeRequest({ method: 'delete', url }),
}

const auth = {
  login: async (form: UserLogin) => await login(form),
  logout: async () => {
    await REQUESTS.post('/auth/logout')
    store.emptyCurrentUser()
  },
  isLoggedIn: async () => {
    // Can't check httpOnly cookie from JavaScript
    // Instead, try to get current user to verify authentication
    const [, error] = await REQUESTS.get('/users/me')
    return error === null
  },
  sendResetPasswordLink: async (email: string) =>
    await REQUESTS.get('/auth/password-forgotten', { user_email: email }),
}

const admin = {
  runWebScraping: async () =>
    await REQUESTS.post('/admin/external-db/web-scraping/run'),
  getWebScrapingState: async () =>
    await REQUESTS.get('/admin/external-db/web-scraping/state'),
  getWebScrapingLastRelease: async () =>
    await REQUESTS.get('/admin/external-db/web-scraping/last-release'),
  getLogFileContent: async () => await REQUESTS.get('/admin/log'),
  getBackupLogFileContent: async (backupNumber: number) =>
    await REQUESTS.get(`/admin/log/backup/${backupNumber}`),
}

const runInfos = {
  list: async (limit: number | null = null, offset: number | null = null) =>
    await REQUESTS.get('/run-infos', { limit, offset }),
  durationStats: async () =>
    await REQUESTS.get('/run-infos/duration-statistics'),
  modeDistribution: async () =>
    await REQUESTS.get('/run-infos/mode-distribution'),
  countPerDay: async () => await REQUESTS.get('/run-infos/count-per-day'),
  sequenceNumberStats: async () =>
    await REQUESTS.get('/run-infos/sequence-number-statistics'),
}

const users = {
  register: async (form: UserForm) => await REQUESTS.post('/users', form),
  me: async () => await REQUESTS.get('/users/me'),
  list: async (limit: number | null = null, offset: number | null = null) =>
    await REQUESTS.get('/users', { limit, offset }),
  updatePassword: async (form: UserPasswordForm, reset: boolean) => {
    return reset
      ? await REQUESTS.put('/users/me/password?reset=true', form)
      : await REQUESTS.put('/users/me/password', form)
  },
  updateProfile: async (form: UserProfileForm) =>
    await REQUESTS.put('/users/me/profile', form),
  updateUserRole: async (id: string, form: UserRoleForm) =>
    await REQUESTS.put(`/users/${id}/role`, form),
  deleteMe: async () => await REQUESTS.delete('/users/me'),
  deleteUser: async (id: string) => await REQUESTS.delete(`/users/${id}`),
}

const codonTables = {
  list: async () => await REQUESTS.get('/codon-tables'),
  get: async (id: string) => await REQUESTS.get(`/users/me/codon-tables/${id}`),
  getTranslations: async (id: string) =>
    await REQUESTS.get(`/users/me/codon-tables/${id}/translations`),
  add: async (form: CodonTableForm) =>
    await REQUESTS.post('/users/me/codon-tables', form),
  update: async (id: string, form: CodonTableForm) =>
    await REQUESTS.put(`/users/me/codon-tables/${id}`, form),
  delete: async (id: string) =>
    await REQUESTS.delete(`/users/me/codon-tables/${id}`),
}

const results = {
  list: async (limit: number | null = null, offset: number | null = null) =>
    await REQUESTS.get(`/users/me/results`, { limit, offset }),
  count: async () => await REQUESTS.get(`/users/me/results/count`),
  get: async (id: string) => await REQUESTS.get(`/users/me/results/${id}`),
  delete: async (id: string) =>
    await REQUESTS.delete(`/users/me/results/${id}`),
  deleteAll: async () => await REQUESTS.delete(`/users/me/results`),
}

const tunedSequences = {
  list: async (resultId: string) =>
    await REQUESTS.get(`/results/${resultId}/tuned-sequences`),
}

export const API = {
  runTraining: async (form: RunTrainingForm) =>
    (await REQUESTS.post('/run-tuning', form)) as ApiResponse<string>,
  cancelTuning: async (jobId: string) =>
    (await REQUESTS.delete(`/tuning/${jobId}`)) as ApiResponse<string>,
  auth: auth,
  users: users,
  codonTables: codonTables,
  results: results,
  tunedSequences: tunedSequences,
  runInfos: runInfos,
  admin: admin,
  compareSequences: async (form: SequenceComparatorForm) =>
    (await REQUESTS.post(
      '/compare-sequences',
      form,
    )) as ApiResponse<TunedSequence>,
}

export async function setCurrentUserInStore() {
  const [data, error] = await API.users.me()

  if (error === null) {
    store.setCurrentUser(data)
  }
}
