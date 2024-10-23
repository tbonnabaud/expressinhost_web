import axios, { type AxiosRequestConfig } from 'axios'
import type { RunTrainingForm, TuningOutput } from './interfaces'

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
      return Promise.reject(error.message)
    } else {
      const message = 'Network error or request failed'
      console.error(message)
      return Promise.reject(message)
    }
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

const REQUESTS = {
  get: async (url: string, params?: object) =>
    await makeRequest({ method: 'get', url, params }),
  post: async (url: string, data?: object) =>
    await makeRequest({ method: 'post', url, data }),
  put: async (url: string, data?: object) =>
    await makeRequest({ method: 'put', url, data }),
  delete: async (url: string) => await makeRequest({ method: 'delete', url }),
}

export const API = {
  runTraining: async (form: RunTrainingForm) =>
    (await REQUESTS.post('/run-tuning', form)) as ApiResponse<TuningOutput>,
}
