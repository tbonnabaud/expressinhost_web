import axios, { type AxiosResponse } from 'axios'
import type { RunTrainingForm, TuningOutput } from './interfaces'

const client = axios.create({
  baseURL: '/api',
})

// client.interceptors.response.use(
//   response => response,
//   error => {
//     console.error(error)
//     return Promise.reject(error)
//   },
// )

function parseResponse(response: AxiosResponse) {
  if (response.status == 200) {
    return response.data
  } else {
    console.error(response.status, response.statusText, response.data)
    return null
  }
}

const REQUESTS = {
  get: async (url: string, params?: object) => {
    const response = await client.get(url, { params })
    return parseResponse(response)
  },
  post: async (url: string, data?: object) => {
    const response = await client.post(url, data)
    return parseResponse(response)
  },
  put: async (url: string, data?: object) => {
    const response = await client.put(url, data)
    return parseResponse(response)
  },
  delete: async (url: string) => {
    const response = await client.delete(url)
    return parseResponse(response)
  },
}

export const API = {
  runTraining: async (form: RunTrainingForm) =>
    (await REQUESTS.post('/run-tuning', form)) as TuningOutput,
}
