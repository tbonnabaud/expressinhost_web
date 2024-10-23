import axios from 'axios'
import type { RunTrainingForm, TuningOutput } from './interfaces'

const client = axios.create({
  baseURL: '/api',
})

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

const REQUESTS = {
  get: async (url: string, params?: object) => {
    try {
      const response = await client.get(url, { params })
      return response.data
    } catch {
      return null
    }
  },
  post: async (url: string, data?: object) => {
    try {
      const response = await client.post(url, data)
      return response.data
    } catch {
      return null
    }
  },
  put: async (url: string, data?: object) => {
    try {
      const response = await client.put(url, data)
      return response.data
    } catch {
      return null
    }
  },
  delete: async (url: string) => {
    try {
      const response = await client.delete(url)
      return response.data
    } catch {
      return null
    }
  },
}

export const API = {
  runTraining: async (form: RunTrainingForm) =>
    (await REQUESTS.post('/run-tuning', form)) as TuningOutput,
}
