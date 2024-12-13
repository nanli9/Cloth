
#include <iostream>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include "shader.h"
#include "camera.h"
#include "Cloth.h"
#include "rigidBody.h"

#define IMGUI_DEFINE_MATH_OPERATORS
#include "imgui.h"
#include "imgui.cpp"
#include "imgui_draw.cpp"
#include "imgui_widgets.cpp"
#include "imgui_tables.cpp"
#include "imgui_demo.cpp"

#include "backends/imgui_impl_opengl3.h"
#include "backends/imgui_impl_opengl3.cpp"
#include "backends/imgui_impl_glfw.h"
#include "backends/imgui_impl_glfw.cpp"

void framebuffer_size_callback(GLFWwindow* window, int width, int height);
void processInput(GLFWwindow* window);
void mouse_callback(GLFWwindow* window, double xpos, double ypos);
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset);
void mouse_button_callback(GLFWwindow* window, int button, int action, int mods);
vec3 unProject(int mouseX, int mouseY, float depth);
mat4 shadowProjection(float a, float b, float c, float d);
// settings
const unsigned int SCR_WIDTH = 800;
const unsigned int SCR_HEIGHT = 600;

// camera
Camera camera(glm::vec3(1.0f, 4.0f, 12.0f));
float lastX = SCR_WIDTH / 2.0f;
float lastY = SCR_HEIGHT / 2.0f;
bool firstMouse = true;

// timing
float deltaTime = 0.0f;
float lastFrame = 0.0f;

// lighting
vec3 lightPos(5.0f, 8.0f, 6.0f);
bool isDragging = false;
bool isDraggingCamera = false;
Cloth c(vec3(0.0, 0 , 0));
vector<RigidBody> rigidBodies;
int main()
{
#ifdef __APPLE__
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

    glfwInit();
    glfwWindowHint(GLFW_SAMPLES, 8);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_COMPAT_PROFILE);
    
    // glfw window creation
    // --------------------
    GLFWwindow* window = glfwCreateWindow(SCR_WIDTH, SCR_HEIGHT, "LearnOpenGL", NULL, NULL);
    if (window == NULL)
    {
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);
    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);

    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    glfwSetCursorPosCallback(window, mouse_callback);
    glfwSetScrollCallback(window, scroll_callback);
    glfwSetMouseButtonCallback(window, mouse_button_callback);
    // tell GLFW to capture our mouse
    //glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);

    // glad: load all OpenGL function pointers
    // ---------------------------------------
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        std::cout << "Failed to initialize GLAD" << std::endl;
        return -1;
    }


    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO();
    ImGui::StyleColorsDark();

    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 330"); 

    // configure global opengl state
    // -----------------------------
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_MULTISAMPLE); // enabled by default on some drivers, but not all so always enable to make sure

    // build and compile our shader 
    // ------------------------------------
    Shader lightingShader("shader/vertex_shader_Phong.glsl", "shader/fragment_shader_Phong.glsl");

    float planeVertices[] = {
        // positions          // texture Coords 
         5.0f, 0.0f,  5.0f,  0.0f, 1.0f,0.0f,
        -5.0f, -0.0f,  5.0f,  0.0f, 1.0f,0.0f,
        -5.0f, -0.0f, -5.0f,  0.0f, 1.0f,0.0f,

         5.0f, -0.0f,  5.0f,  0.0f, 1.0f,0.0f,
        -5.0f, -0.0f, -5.0f,  0.0f, 1.0f,0.0f,
         5.0f, -0.0f, -5.0f,  0.0f, 1.0f,0.0f,
    };

    //create VAO
    unsigned int planeVAO, planeVBO;
    glGenVertexArrays(1, &planeVAO);
    glGenBuffers(1, &planeVBO);
    glBindVertexArray(planeVAO);
    glBindBuffer(GL_ARRAY_BUFFER, planeVBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(planeVertices), &planeVertices, GL_STATIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(3 * sizeof(float)));
    RigidBody r;
    rigidBodies.push_back(r);

    lightingShader.use();
    // render loop
    // -----------
    while (!glfwWindowShouldClose(window))
    {
        float currentFrame = static_cast<float>(glfwGetTime());
        deltaTime = currentFrame - lastFrame;
        lastFrame = currentFrame;

        // input
        // -----
        processInput(window);
        
       /* double xpos, ypos;
        glfwGetCursorPos(window, &xpos, &ypos);
        std::cout << "Cursor position: " << xpos << ", " << ypos << std::endl;*/
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        ImGui::Begin("Cloth Parameters");
        ImGui::SliderFloat("Stretch", &c.stretchCompliance, 0, 0.1f);
        ImGui::SliderFloat("Bend", &c.bendCompliance, 0, 0.1f);
        ImGui::SliderFloat("K_damping", &c.k_damping, 0.0, 1.0f);
        if (ImGui::Button("Wind")) {
            c.f_external.x = -5.0f;
            c.f_external.z = -3.0f;
        }
        if (ImGui::Button("gravity")) {
            c.f_external.y = -9.8f;
        }
        if (ImGui::Button("reset external forces")) {
            c.f_external = vec3(0,0,0);
        }
        if (ImGui::Button("reset")) {
            c.reset();
        }
        if (ImGui::Button("pin")) {
            c.pinned = !c.pinned;
        }
        if (ImGui::Button("Line Mode")) {
            c.lineDisplay = !c.lineDisplay;
        }
        if (ImGui::Button("self-collision")) {
            c.handleSelfCollision = !c.handleSelfCollision;
        }
        ImGui::End();

        // Rendering
        ImGui::Render();

        //glClearColor(0.529f, 0.808f, 0.922f, 1.0f);
        glClearColor(0, 0, 0, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        mat4 projection = glm::perspective(glm::radians(camera.Zoom), (float)SCR_WIDTH / (float)SCR_HEIGHT, 0.1f, 100.0f);
        mat4 view = camera.GetViewMatrix();
        mat4 model = mat4(1.0f);
        model = translate(model, vec3(0,-30,0));
        model = scale(model, vec3(20, 20, 20));
        lightingShader.setMat4("projection", projection);
        lightingShader.setMat4("view", view);

        lightingShader.setVec3("lightPos", lightPos);
        lightingShader.setVec3("eye", camera.Position);
        lightingShader.setVec3("color", vec3(0.827f, 0.827f, 0.827f));

        glBindVertexArray(planeVAO);
        lightingShader.setMat4("model", model);
        glDrawArrays(GL_TRIANGLES, 0, 6);
        glBindVertexArray(0);

        mat4 shadowMatrix = shadowProjection(0, 1, 0, 2.9);

        model = mat4(1.0f);
        /*model = translate(model, vec3(2, -1, 2));
        model = scale(model, vec3(1, 1, 1));*/
        lightingShader.setMat4("model", model);
        lightingShader.setVec3("color", vec3(0.0f, 1.0f, 0.0f));

        for(auto& r: rigidBodies)
            r.draw(lightingShader, shadowMatrix);

        
        model = mat4(1.0f);
        model = translate(model, vec3(0, 0, 0));
        lightingShader.setMat4("model", model);

        c.draw(lightingShader, shadowMatrix);
        c.update(std::min(deltaTime,0.03f), rigidBodies);

        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

        // glfw: swap buffers and poll IO events (keys pressed/released, mouse moved etc.)
        // -------------------------------------------------------------------------------
        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    // glfw: terminate, clearing all previously allocated GLFW resources.
    // ------------------------------------------------------------------
    glfwTerminate();

    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
    glfwDestroyWindow(window);
    glfwTerminate();

    return 0;
}

mat4 shadowProjection(float a, float b, float c, float d)
{
    mat4 M;
    vec4 L = vec4(lightPos, 1.0);
    float lambda = dot(L, vec4(a,b,c,d));
    //first row
    M[0][0] = lambda - a * lightPos.x;
    M[1][0] = -a * lightPos.y;
    M[2][0] = -a * lightPos.z;
    M[3][0] = -a;
    //second row
    M[0][1] = -b * lightPos.x;
    M[1][1] = lambda - b * lightPos.y;
    M[2][1] = -b * lightPos.z;
    M[3][1] = -b;
    //third row
    M[0][2] = -c * lightPos.x;
    M[1][2] = -c * lightPos.y;
    M[2][2] = lambda - c * lightPos.z;
    M[3][2] = -c;
    //fourth row
    M[0][3] = -d * lightPos.x;
    M[1][3] = -d * lightPos.y;
    M[2][3] = -d * lightPos.z;
    M[3][3] = lambda - d;

    return transpose(M);
}


// process all input: query GLFW whether relevant keys are pressed/released this frame and react accordingly
// ---------------------------------------------------------------------------------------------------------
void processInput(GLFWwindow* window)
{
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);

    if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
        camera.ProcessKeyboard(FORWARD, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
        camera.ProcessKeyboard(BACKWARD, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
        camera.ProcessKeyboard(LEFT, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
        camera.ProcessKeyboard(RIGHT, deltaTime);

    if (glfwGetKey(window, GLFW_KEY_LEFT) == GLFW_PRESS)
        camera.ProcessKeyboardRotation(ARROW_LEFT, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_RIGHT) == GLFW_PRESS)
        camera.ProcessKeyboardRotation(ARROW_RIGHT, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_UP) == GLFW_PRESS)
        camera.ProcessKeyboardRotation(ARROW_UP, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_DOWN) == GLFW_PRESS)
        camera.ProcessKeyboardRotation(ARROW_DOWN, deltaTime);

    if (glfwGetKey(window, GLFW_KEY_R) == GLFW_PRESS)
        c.f_external.x = 0.0f;
}

// glfw: whenever the window size changed (by OS or user resize) this callback function executes
// ---------------------------------------------------------------------------------------------
void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    // make sure the viewport matches the new window dimensions; note that width and 
    // height will be significantly larger than specified on retina displays.
    glViewport(0, 0, width, height);
}

void mouse_callback(GLFWwindow* window, double xposIn, double yposIn)
{
    float xpos = static_cast<float>(xposIn);
    float ypos = static_cast<float>(yposIn);

    if (firstMouse)
    {
        lastX = xpos;
        lastY = ypos;
        firstMouse = false;
    }

    float xoffset = xpos - lastX;
    float yoffset = lastY - ypos; // reversed since y-coordinates go from bottom to top
    lastX = xpos;
    lastY = ypos;
    if (isDragging)
    {
        
        double mouseX, mouseY;
        glfwGetCursorPos(window, &mouseX, &mouseY);
        float depth;
        glReadPixels(static_cast<int>(mouseX), static_cast<int>(SCR_HEIGHT - mouseY), 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &depth);
        double x, y, z;
        vec3 Pos = unProject(static_cast<int>(mouseX), static_cast<int>(mouseY), depth);
        //c.f_external.x = std::min(xoffset,2.0f);
        if (c.grabIndex != -1)
        {
            c.particles[c.grabIndex].x += 0.01f * vec3(xoffset, yoffset,0);
            c.particles[c.grabIndex].p += 0.01f * vec3(xoffset, yoffset, 0);
        }
        else
            c.grab(Pos);
    }
    if (isDraggingCamera)
        camera.ProcessMouseMovement(xoffset, yoffset);
    
    
}
void mouse_button_callback(GLFWwindow* window, int button, int action, int mods) {
    if (action == GLFW_PRESS) {
        if (button == GLFW_MOUSE_BUTTON_LEFT) {
            isDragging = true;
        }
        else if (button == GLFW_MOUSE_BUTTON_RIGHT) {
            /*double mouseX, mouseY;
            glfwGetCursorPos(window, &mouseX, &mouseY);
            std::cout << "Right mouse button pressed at (" << mouseX << ", " << mouseY << ")\n";*/
            isDraggingCamera = true;
        }
    }
    else if (action == GLFW_RELEASE)
    {
        if (button == GLFW_MOUSE_BUTTON_LEFT) {
            isDragging = false;
            if (c.grabIndex != -1)
            {
                c.particles[c.grabIndex].inverseMass = c.grabPointInvereMass;
                c.grabIndex = -1;
            }
        }
        else if (button == GLFW_MOUSE_BUTTON_RIGHT) {
            isDraggingCamera = false;
        }
        
    }
}
// glfw: whenever the mouse scroll wheel scrolls, this callback is called
// ----------------------------------------------------------------------
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
{
    camera.ProcessMouseScroll(static_cast<float>(yoffset));
}
vec3 unProject(int mouseX, int mouseY, float depth)
{
    vec3 pos;
    mat4 projection = glm::perspective(glm::radians(camera.Zoom), (float)SCR_WIDTH / (float)SCR_HEIGHT, 0.1f, 100.0f);
    mat4 view = camera.GetViewMatrix();
    mat4 model = mat4(1.0f);
    float x = (2.0f * mouseX) / SCR_WIDTH - 1.0f;
    float y = 1.0f - (2.0f *  mouseY) / SCR_HEIGHT;
    vec4 clipCoord(x, y, 2*depth-1,1.0f);
    vec4 projectCoord =  inverse(projection) * clipCoord;
    vec4 eyeCoord = inverse(view) * projectCoord;
    vec4 worldCoord = inverse(model) * eyeCoord;

    return vec3(worldCoord) / worldCoord.w;
}