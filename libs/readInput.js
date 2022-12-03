const fs = require('fs/promises')
const path = require('path')


const readInput = async (dir, inputFile, splitter = /\r?\n/) => {
    const data = await fs.readFile(path.resolve(dir, inputFile), 'utf8')
    return data.trim().split(splitter)
}

module.exports = readInput
