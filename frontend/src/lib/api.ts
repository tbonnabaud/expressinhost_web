import axios, { type AxiosRequestConfig } from 'axios'
import type { RunTrainingForm, TuningOutput } from './interfaces'

const client = axios.create({
  baseURL: '/api',
})

// Interceptor to unify error handling
client.interceptors.response.use(
  response => response,
  error => {
    if (axios.isAxiosError(error) && error.response) {
      console.error(error.message, error.response.data)
    } else {
      // 'Network error or request failed'
      console.error(error.message)
    }

    return Promise.reject(error)
  },
)

async function request(config: AxiosRequestConfig) {
  try {
    const response = await client.request(config)
    return response.data
  } catch {
    return null
  }
}

const REQUESTS = {
  get: async (url: string, params?: object) =>
    request({ method: 'get', url: url, params: params }),
  post: async (url: string, data?: object) =>
    request({ method: 'post', url: url, data: data }),
  put: async (url: string, data?: object) =>
    request({ method: 'put', url: url, data: data }),
  delete: async (url: string) => request({ method: 'delete', url: url }),
}

export const API = {
  runTraining: async (form: RunTrainingForm) =>
    (await REQUESTS.post('/run-tuning', form)) as TuningOutput,
}
