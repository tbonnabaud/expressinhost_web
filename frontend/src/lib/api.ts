import axios, { type AxiosRequestConfig } from 'axios'
import type {
  RunTrainingForm,
  TuningOutput,
  UserForm,
  UserLogin,
  Token,
} from './interfaces'

type ApiResponse<T> = [T | null, string | null]

const client = axios.create({
  baseURL: '/api',
})

// Interceptor to unify error handling
client.interceptors.response.use(
  response => response,
  error => {
    if (axios.isAxiosError(error) && error.response) {
      console.error(error.message, error.response.data)

      if (error.response.status == 401) {
        localStorage.getItem('accessToken')
      }

      return Promise.reject(error.message)
    } else {
      const message = 'Network error or request failed'
      console.error(message)
      return Promise.reject(message)
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
    return [null, error]
  }
}

async function postForm(url: string, form: UserLogin) {
  try {
    const response = await client.postForm(url, form)
    return [response.data, null] as [Token, null]
  } catch (error) {
    return [null, error] as [null, string]
  }
}

async function login(form: UserLogin) {
  const [data, error] = await postForm('/token', form)

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

const users = {
  register: async (form: UserForm) => await REQUESTS.post('/users', form),
  login: async (form: UserLogin) => await login(form),
  logout: () => localStorage.removeItem('accessToken'),
  isLoggedIn: () => localStorage.getItem('accessToken') !== null,
  me: async () => await REQUESTS.get('/users/me'),
}

export const API = {
  runTraining: async (form: RunTrainingForm) =>
    (await REQUESTS.post('/run-tuning', form)) as ApiResponse<TuningOutput>,
  users: users,
}
