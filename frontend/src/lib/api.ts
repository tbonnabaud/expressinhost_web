import axios, { type AxiosResponse } from 'axios'
import type { RunTrainingForm, TuningOutput } from './interfaces'

const client = axios.create({
  baseURL: '/api',
})

function parseResponse(response: AxiosResponse) {
  if (response.status == 200) {
    return response.data
  } else {
    console.error(response.status, response.statusText, response.data)
    return null
  }
}

export const API = {
  runTraining: async (form: RunTrainingForm) => {
    const response = await client.post('/run-tuning', form)
    return parseResponse(response) as TuningOutput
  },
}
