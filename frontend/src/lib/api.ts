import axios, { type AxiosRequestConfig } from 'axios'
import type {
  RunTrainingForm,
  TuningOutput,
  UserForm,
  UserLogin,
  Token,
  CodonTableForm,
  UserPasswordForm,
  UserProfileForm,
  UserRoleForm,
} from './interfaces'
import { store } from './store'

interface ApiError {
  code: number | null
  detail: string
}

type ApiResponse<T> = [T, null] | [null, ApiError]

const client = axios.create({
  baseURL: '/api',
})

// Interceptor to unify error handling
client.interceptors.response.use(
  response => response,
  error => {
    if (axios.isAxiosError(error) && error.response) {
      if (error.response.status == 401) {
        console.warn(error.message, error.response.data)
        localStorage.removeItem('accessToken')
        store.emptyCurrentUser()
      } else {
        console.error(error.message, error.response.data)
      }

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

// Interceptor to add access token if exists
client.interceptors.request.use(
  request => {
    const accessToken = localStorage.getItem('accessToken')
    if (accessToken) {
      request.headers['Authorization'] = `Bearer ${accessToken}`
    }
    return request
  },
  error => {
    console.error(error)
    return Promise.reject(error)
  },
)

async function makeRequest(config: AxiosRequestConfig) {
  try {
    const response = await client.request(config)
    return [response.data, null]
  } catch (error) {
    return [null, error] as [null, ApiError]
  }
}

async function postForm(url: string, form: UserLogin) {
  try {
    const response = await client.postForm(url, form)
    return [response.data, null] as [Token, null]
  } catch (error) {
    return [null, error] as [null, ApiError]
  }
}

async function login(form: UserLogin) {
  const [data, error] = await postForm('/auth/token', form)

  if (!error && data) {
    localStorage.setItem('accessToken', data.access_token)
  }

  return error
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
  logout: () => localStorage.removeItem('accessToken'),
  isLoggedIn: () => localStorage.getItem('accessToken') !== null,
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
    (await REQUESTS.post('/run-tuning', form)) as ApiResponse<TuningOutput>,
  auth: auth,
  users: users,
  codonTables: codonTables,
  results: results,
  tunedSequences: tunedSequences,
  runInfos: runInfos,
  admin: admin,
}

export async function setCurrentUserInStore() {
  const [data, error] = await API.users.me()

  if (error === null) {
    store.setCurrentUser(data)
  }
}
