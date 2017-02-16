    viz->vertex_array_buffer.set_buffer_data(3 * sizeof(GLfloat)
        * m_vertices.size(), m_vertices.front().data());
    viz->normal_array_buffer.set_buffer_data(3 * sizeof(GLfloat)
        * m_normals.size(), m_normals.front().data());
    viz->index_array_buffer.set_buffer_data(3 * sizeof(GLuint)
        * m_faces.size(), m_faces.front().data());
