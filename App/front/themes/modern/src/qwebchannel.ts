/* eslint-disable @typescript-eslint/no-explicit-any */
export class QWebChannel {
    transport: any;
    execCallbacks: { [key: number]: (data: any) => void };
    execId: number;
    objects: { [key: string]: any };

    constructor(transport: any, initCallback: (channel: QWebChannel) => void) {
        this.transport = transport;
        this.execCallbacks = {};
        this.execId = 0;
        this.objects = {};

        this.transport.onmessage = (message: any) => {
            const data = JSON.parse(message.data);

            if (data.type === 1) { // Properties update
                for (const i in data.properties) {
                    const prop = data.properties[i];
                    if (this.objects[prop.object]) {
                        this.objects[prop.object].wrappedProperties[prop.property] = prop.value;
                    }
                }
            } else if (data.type === 2) { // Method response
                const callback = this.execCallbacks[data.id];
                if (callback) {
                    delete this.execCallbacks[data.id];
                    callback(data.response);
                }
            } else if (data.type === 5) { // Initialization payload containing backend objects mappings
                for (const name in data.objectData) {
                    this.objects[name] = this.createQtObject(name, data.objectData[name]);
                }
                // Invoke callback only after objects are fully registered and built
                initCallback(this);
            }
        };

        this.transport.send(JSON.stringify({ type: 0 })); // Request initialization handshake
    }

    createQtObject(name: string, signalMethodMap: any) {
        const obj: any = { wrappedProperties: {} };
        // Map methods directly to invocation packets sent down the socket transport
        if (signalMethodMap.methods) {
            for (let i = 0; i < signalMethodMap.methods.length; i++) {
                const methodName = signalMethodMap.methods[i][0];
                obj[methodName] = (...args: any[]) => {
                    let callback: ((data: any) => void) | null = null;
                    if (args.length > 0 && typeof args[args.length - 1] === 'function') {
                        callback = args.pop();
                    }
                    this.execId++;
                    if (callback) {
                        this.execCallbacks[this.execId] = callback;
                    }
                    this.transport.send(JSON.stringify({
                        type: 3,
                        id: this.execId,
                        object: name,
                        method: methodName,
                        args: args
                    }));
                };
            }
        }
        return obj;
    }
}